#' GFISHER reef-habitat basemaps for WFS Ecospace.
#'
#' Builds per-cell proportional coverage of six reef classes from the FWRI
#' East Gulf side-scan geodatabase:
#'
#'   AL / AM / AH  artificial reef, low / medium / high relief
#'   NL / NM / NH  natural reef,    low / medium / high relief
#'
#' Coverage is only *observed* where side-scan mapping exists (the microgrid
#' footprint). Elsewhere the natural classes are extrapolated and the
#' artificial classes are set to zero -- artificial structure is placed, not
#' gradient-distributed, so predicting it into unsurveyed water is not
#' defensible.
#'
#' Everything here is terra + sf. Do not add library(raster) or library(sp) to
#' this file: the driver attaches terra first and sources R/ afterwards, so a
#' library() call here would mask terra's generics for the whole session.
#'
#' Sampling-bias caveat: side-scan effort is targeted rather than random, so
#' mapped cells are probably reef-enriched relative to the shelf as a whole.
#' Every fill method inherits that bias; none of them can correct it.


#' Fill unmapped cells in a habitat proportion raster.
#'
#' @param r SpatRaster of observed proportion; NA outside the mapped footprint.
#' @param depth SpatRaster of positive depths, the template grid.
#' @param mapped Logical SpatRaster, TRUE where side-scan mapping exists.
#' @param max.depth Only fill cells at or shallower than this. Inf fills all.
#' @param method One of "gam", "strata", "idw", "krige", "none".
#' @param covariates Optional SpatRaster stack on the depth grid (e.g. the
#'   dbSEABED substrate fractions). Used by "gam" only.
#' @param strata SpatRaster of stratum IDs on the depth grid, used by "strata".
#'   NULL bins `depth` at `strata.breaks`, the only stratification available for
#'   cells the survey never visited.
#' @param strata.breaks Depth breaks (m) used when `strata` is NULL.
#' @param idw.idp,idw.nmax Inverse-distance power and neighbour count.
#' @param k.xy Basis dimension for the spatial smooth in the GAM.
#' @param max.cov.na Drop a covariate whose NA fraction exceeds this.
#' @param anchor.zero Add synthetic z = 0 observations to the training set:
#'   "land" (cells where depth is NA), "deep" (cells past `max.depth`), "both",
#'   or "none" (default). This reproduces the legacy behaviour, where zeroing
#'   land before the fill meant those cells entered gstat as real observations.
#'   They are not observations -- nobody surveyed them -- and because effort is
#'   concentrated inshore where reef is densest, including them biases the fill
#'   downward there (NL filled mean 0.0118 -> 0.0049 at 5 min). Distance-based
#'   methods only; ignored with a warning for "gam" and "strata".
#' @param clamp.to.observed Cap filled values at the largest observed value for
#'   this class. Stops a model extrapolating past anything the survey ever saw.
#' @param verbose Print fit diagnostics.
#' @return SpatRaster with filled values at eligible cells; other cells unchanged.
fn.fill_habitat_gaps <- function(r, depth, mapped, max.depth = Inf,
                                 method = c("gam", "strata", "idw", "krige", "none"),
                                 covariates = NULL,
                                 strata = NULL,
                                 strata.breaks = c(0, 20, 40, 60, 100, 200, Inf),
                                 idw.idp = 4, idw.nmax = 8,
                                 k.xy = 60, max.cov.na = 0.5,
                                 anchor.zero = c("none", "land", "deep", "both"),
                                 clamp.to.observed = TRUE, verbose = TRUE) {


  # r=prop.i
  # max.depth = max.depth.hab
  # method = fill.method
  # k.xy = 60
  # max.cov.na = 0.5

  method <- match.arg(method,c("gam", "strata", "idw", "krige", "none"))
  anchor.zero <- match.arg(anchor.zero, c("none", "land", "deep", "both"))
  if (method == "none") return(r)

  # Equal-area CRS for anything distance-based. Lon/lat degrees are not a
  # sensible metric for IDW or a variogram at this latitude.
  ea.crs <- paste("+proj=aea +lat_1=24 +lat_2=31.5 +lat_0=23 +lon_0=-84",
                  "+datum=WGS84 +units=m +no_defs")

  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  d  <- data.frame(cell   = seq_len(terra::ncell(r)),
                   x      = xy[, 1],
                   y      = xy[, 2],
                   z      = terra::values(r, mat = FALSE),
                   dep    = terra::values(depth, mat = FALSE),
                   mapped = as.logical(terra::values(mapped, mat = FALSE)),
                   wt     = 1)

  # Weight observations by how much of the cell was actually surveyed: a cell
  # with two microgrids scanned says far less than one with two thousand.
  if (!is.null(attr(r, "mapped.area"))) d$wt <- attr(r, "mapped.area")

  cov.names <- character(0)
  if (!is.null(covariates)) {
    if (!terra::compareGeom(covariates, depth, stopOnError = FALSE)) {
      warning("covariates do not share the depth geometry; resampling.")
      covariates <- terra::resample(covariates, depth, method = "bilinear")
    }
    cv <- as.data.frame(terra::values(covariates, mat = TRUE))
    names(cv) <- make.names(names(covariates), unique = TRUE)
    cov.names <- names(cv)
    d <- cbind(d, cv)
  }

  train <- d[d$mapped & !is.na(d$z), , drop = FALSE]
  pred  <- d[!d$mapped & !is.na(d$dep) & d$dep <= max.depth, , drop = FALSE]

  # Optionally re-create the legacy behaviour, where setting land (and, in the
  # nested variant, everything past the depth cut-off) to 0 before the fill
  # meant those cells entered gstat as genuine z = 0 observations. They are not
  # observations -- nobody surveyed them -- so they are excluded by default.
  # Distance-based methods only: land rows have dep = NA and would be dropped
  # by the GAM anyway.
  if (anchor.zero != "none") {
    if (!method %in% c("idw", "krige")) {
      warning("anchor.zero is ignored for method = '", method, "'.")
    } else {
      is.land <- is.na(d$dep)
      is.deep <- !is.na(d$dep) & d$dep > max.depth
      keep <- switch(anchor.zero,
                     land = is.land,
                     deep = is.deep,
                     both = is.land | is.deep)
      add <- d[keep, , drop = FALSE]
      if (nrow(add) > 0) {
        add$z  <- 0
        add$wt <- 1
        train  <- rbind(train, add)
        if (verbose)
          message(sprintf("    anchor.zero='%s': added %d synthetic zero(s) to %d mapped cells",
                          anchor.zero, nrow(add), sum(d$mapped & !is.na(d$z))))
      }
    }
  }

  if (nrow(pred) == 0L) {
    if (verbose) message("    nothing to fill.")
    return(r)
  }
  if (nrow(train) < 30L) {
    warning("Only ", nrow(train), " mapped cells; skipping fill.")
    return(r)
  }

  fit.vals <- switch(
    method,

    gam = {
      if (!requireNamespace("mgcv", quietly = TRUE))
        stop("fill.method = 'gam' needs the mgcv package.")

      # Keep covariates with decent coverage and enough spread to smooth over.
      # dbSEABED has scattered gaps, so a strict no-NA rule would reject every
      # substrate layer; instead drop the badly incomplete ones and impute the
      # rest at the training median so no prediction row is lost.
      use.cov <- character(0)
      for (nm in cov.names) {
        frac.na <- mean(is.na(c(train[[nm]], pred[[nm]])))
        if (frac.na > max.cov.na) {
          if (verbose) message(sprintf("    skipping covariate %s (%.0f%% NA)",
                                       nm, 100 * frac.na))
          next
        }
        if (length(unique(stats::na.omit(train[[nm]]))) < 10L) next
        med <- stats::median(train[[nm]], na.rm = TRUE)
        train[[nm]][is.na(train[[nm]])] <- med
        pred[[nm]][is.na(pred[[nm]])]   <- med
        use.cov <- c(use.cov, nm)
      }

      k.xy <- min(k.xy, max(10L, floor(nrow(train) / 4)))
      terms <- c(sprintf("s(x, y, k = %d)", k.xy), "s(dep, k = 5)",
                 sprintf("s(%s, k = 5)", use.cov))
      form <- stats::as.formula(paste("z ~", paste(terms, collapse = " + ")))

      # Quasibinomial with area weights: the response is a proportion of area,
      # bounded [0, 1] with many exact zeros. Beta cannot take 0 or 1, and
      # Gaussian/Tweedie ignore the bound. Weights carry survey effort.
      w <- train$wt / mean(train$wt, na.rm = TRUE)
      fit <- try(mgcv::gam(form, data = train, family = stats::quasibinomial(),
                           weights = w, method = "REML"), silent = TRUE)

      if (inherits(fit, "try-error")) {
        warning("GAM failed (", conditionMessage(attr(fit, "condition")),
                "); falling back to IDW.")
        return(fn.fill_habitat_gaps(r, depth = depth, mapped = mapped,
                                    max.depth = max.depth, method = "idw",
                                    covariates = covariates, idw.idp = idw.idp,
                                    idw.nmax = idw.nmax, k.xy = k.xy,
                                    max.cov.na = max.cov.na, verbose = verbose))
      }
      if (verbose) {
        s <- summary(fit)
        message(sprintf("    gam: n=%d  dev.expl=%.1f%%  terms: %s",
                        nrow(train), 100 * s$dev.expl,
                        paste(c("s(x,y)", "s(dep)", use.cov), collapse = ", ")))
      }
      as.numeric(stats::predict(fit, newdata = pred, type = "response"))
    },

    strata = {
      # Design-based ratio estimator: within each stratum, total habitat area
      # over total scanned area. That is the scanned-area-weighted mean of the
      # per-cell proportions, and it is what the survey's own stratification
      # supports. Unmapped cells in a stratum all take its ratio, so the result
      # is flat within strata by construction.
      s <- strata
      if (is.null(s)) {
        rcl <- cbind(utils::head(strata.breaks, -1), strata.breaks[-1],
                     seq_len(length(strata.breaks) - 1))
        s <- terra::classify(depth, rcl, include.lowest = TRUE, right = TRUE)
      }
      sv <- terra::values(s, mat = FALSE)
      train$stratum <- sv[train$cell]
      pred$stratum  <- sv[pred$cell]

      num <- tapply(train$z * train$wt, train$stratum, sum, na.rm = TRUE)
      den <- tapply(train$wt,           train$stratum, sum, na.rm = TRUE)
      mu  <- num / den
      mu.all <- sum(train$z * train$wt, na.rm = TRUE) / sum(train$wt, na.rm = TRUE)

      if (verbose) {
        tab <- data.frame(stratum = names(mu),
                          n_mapped = as.integer(table(train$stratum)[names(mu)]),
                          ratio    = signif(as.numeric(mu), 3))
        message("    strata means (scanned-area weighted):")
        message(paste(utils::capture.output(print(tab, row.names = FALSE)),
                      collapse = "\n"))
      }

      v <- as.numeric(mu[as.character(pred$stratum)])
      n.miss <- sum(is.na(v))
      if (n.miss > 0 && verbose)
        message(sprintf("    %d cell(s) in strata with no mapped data; using the overall ratio %.5f",
                        n.miss, mu.all))
      v[is.na(v)] <- mu.all
      v
    },

    idw = {
      if (!requireNamespace("gstat", quietly = TRUE))
        stop("fill.method = 'idw' needs the gstat package.")
      tr <- sf::st_transform(sf::st_as_sf(train, coords = c("x", "y"), crs = 4326), ea.crs)
      pr <- sf::st_transform(sf::st_as_sf(pred,  coords = c("x", "y"), crs = 4326), ea.crs)
      g  <- gstat::idw(z ~ 1, locations = tr, newdata = pr,
                       idp = idw.idp, nmax = idw.nmax, debug.level = 0)
      if (verbose) message(sprintf("    idw: n=%d  idp=%g  nmax=%d",
                                   nrow(train), idw.idp, idw.nmax))
      g$var1.pred
    },

    krige = {
      if (!requireNamespace("gstat", quietly = TRUE))
        stop("fill.method = 'krige' needs the gstat package.")
      tr <- sf::st_transform(sf::st_as_sf(train, coords = c("x", "y"), crs = 4326), ea.crs)
      pr <- sf::st_transform(sf::st_as_sf(pred,  coords = c("x", "y"), crs = 4326), ea.crs)
      v  <- gstat::variogram(z ~ 1, tr)
      vm <- try(gstat::fit.variogram(v, gstat::vgm(c("Exp", "Sph", "Gau"))),
                silent = TRUE)
      if (inherits(vm, "try-error")) {
        warning("Variogram fit failed; falling back to IDW.")
        return(fn.fill_habitat_gaps(r, depth = depth, mapped = mapped,
                                    max.depth = max.depth, method = "idw",
                                    covariates = covariates, idw.idp = idw.idp,
                                    idw.nmax = idw.nmax, k.xy = k.xy,
                                    max.cov.na = max.cov.na, verbose = verbose))
      }
      if (verbose) message(sprintf("    krige: n=%d  model=%s  range=%.0f m",
                                   nrow(train), as.character(vm$model[nrow(vm)]),
                                   vm$range[nrow(vm)]))
      gstat::krige(z ~ 1, tr, pr, model = vm, debug.level = 0)$var1.pred
    }
  )

  fit.vals <- as.numeric(fit.vals)

  # A logit-scale smooth can predict well past its training range where the
  # footprint is sparse. Nothing downstream should see more reef in a cell than
  # the survey ever recorded in one.
  if (isTRUE(clamp.to.observed)) {
    obs.max <- max(train$z, na.rm = TRUE)
    n.cap <- sum(fit.vals > obs.max, na.rm = TRUE)
    if (n.cap > 0 && verbose)
      message(sprintf("    clamping %d filled cell(s) to the observed max %.4f",
                      n.cap, obs.max))
    fit.vals <- pmin(fit.vals, obs.max)
  }

  out <- r
  out[pred$cell] <- pmin(pmax(fit.vals, 0), 1)
  out
  #plot(out)
}


#' Build GFISHER reef-habitat proportion maps on the Ecospace grid.
#'
#' @param depth SpatRaster of positive depths (the template grid).
#' @param file.gdb Path to the FWRI East Gulf geodatabase.
#' @param dir.maps Output directory. ASCII and PNGs are written here directly.
#' @param hab.classes Habitat classes to build, in output order.
#' @param drop.classes Values of NewHabStrat to discard outright (e.g. "AP").
#' @param max.depth.hab Deepest cell eligible for extrapolation, in metres.
#'   Pass the model's exclusion depth. NULL means no limit.
#' @param extrapolate Classes to gap-fill. Others get 0 outside the footprint.
#' @param fill.method Gap-fill method; see fn.fill_habitat_gaps().
#' @param covariates Optional SpatRaster stack on the depth grid, used by the
#'   GAM fill (e.g. the dbSEABED substrate fractions).
#' @param strata,strata.breaks Passed to the strata fill. `strata` is a
#'   SpatRaster of stratum IDs; NULL bins depth at `strata.breaks`.
#' @param idw.idp,idw.nmax Passed to the IDW fill.
#' @param anchor.zero Passed to the fill. "none" (default) trains only on cells
#'   the survey actually scanned. "land", "deep" or "both" add synthetic zeros,
#'   reproducing the legacy behaviour; see fn.fill_habitat_gaps() for why that
#'   biases the fill downward inshore. Use
#'   `anchor.zero = "both", max.depth.hab = 200` to match the old maps.
#' @param clamp.to.observed Cap filled values at the largest observed value for
#'   the class.
#' @param cap.total Rescale cells where the six classes sum above 1, so the
#'   downstream sum-to-1 step keeps headroom for seagrass and seabed.
#' @param run.checks Run the microgrid/habitat join QC. Off by default: the
#'   join is used by nothing else and costs a 186k x 309k spatial intersect.
#' @param write.ascii Write Ecospace ASCII grids.
#' @param do.plot Draw progress plots to the active device.
#' @return list(habitat = SpatRaster of the classes, microgrid = mapped area m2).
fn.make_GFISHER_habitat_maps <- function(depth, file.gdb, dir.maps,
                                         hab.classes   = c('AL','AM','AH','NL','NM','NH'),
                                         drop.classes  = 'AP',
                                         max.depth.hab = NULL,
                                         extrapolate   = c('NL','NM','NH'),
                                         fill.method   = c('gam','strata','idw','krige','none'),
                                         covariates    = NULL,
                                         strata        = NULL,
                                         strata.breaks = c(0, 20, 40, 60, 100, 200, Inf),
                                         idw.idp = 4, idw.nmax = 8,
                                         anchor.zero = c("none", "land", "deep", "both"),
                                         clamp.to.observed = TRUE,
                                         cap.total   = TRUE,
                                         run.checks  = FALSE,
                                         write.ascii = TRUE,
                                         do.plot     = TRUE) {
  
  # file.gdb      = file.gfishergdb
  # dir.maps      = file.path(dir.habitats,'gfisher')
  # max.depth.hab = 200
  # hab.classes   = c('AL','AM','AH','NL','NM','NH')
  # drop.classes  = 'AP'
  # extrapolate   = c('NL','NM','NH')
  # fill.method   = c('gam','strata','idw','krige','none')[3]
  # covariates    = NULL
  # strata        = NULL
  # strata.breaks = c(0, 20, 40, 60, 100, 200, Inf)
  # idw.idp = 4
  # idw.nmax = 8
  # clamp.to.observed = TRUE
  # cap.total   = TRUE
  # run.checks  = FALSE
  # write.ascii = FALSE
  # do.plot     = TRUE
  # anchor.zero = 'both'

  fill.method <- match.arg(fill.method,c('gam','strata','idw','krige','none'))

  # ---- template grid --------------------------------------------------------
  if (inherits(depth, "Raster")) depth <- terra::rast(depth)
  stopifnot(inherits(depth, "SpatRaster"))
  # marmap-derived grids often arrive with no CRS; downstream projection needs one.
  if (is.na(terra::crs(depth)) || !nzchar(terra::crs(depth)))
    terra::crs(depth) <- "EPSG:4326"

  res.min <- round(terra::res(depth)[1] * 60, 0)
  if (is.null(max.depth.hab)) max.depth.hab <- Inf
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)

  # ---- read the geodatabase -------------------------------------------------
  # Two layers matter:
  #   *Microgrid_Mapped* - 0.1 x 0.1 nm cells actually scanned (the footprint)
  #   *Hab_Data_Dissolve* - digitised habitat polygons, classed by NewHabStrat
  lyrs <- sf::st_layers(file.gdb)
  nm.micro <- grep("Microgrid", lyrs$name, ignore.case = TRUE, value = TRUE)
  nm.hab   <- grep("Dissolve",  lyrs$name, ignore.case = TRUE, value = TRUE)
  if (length(nm.micro) != 1L)
    stop("Expected one Microgrid layer in the gdb, found ", length(nm.micro), ".")
  if (length(nm.hab) != 1L)
    stop("Expected one Dissolve layer in the gdb, found ", length(nm.hab), ".")

  message("Reading geodatabase...")
  microgrid <- sf::st_read(file.gdb, layer = nm.micro, quiet = TRUE)
  habitat   <- sf::st_read(file.gdb, layer = nm.hab,   quiet = TRUE)
  message(sprintf("  microgrid: %s features   habitat: %s features",
                  format(nrow(microgrid), big.mark = ","),
                  format(nrow(habitat),   big.mark = ",")))

  # ---- habitat classes ------------------------------------------------------
  if (length(drop.classes))
    habitat <- habitat[!habitat$NewHabStrat %in% drop.classes, ]

  # A stray MULTISURFACE record has appeared in past editions and has no
  # centroid method; drop anything that is not a polygon.
  gt  <- as.character(sf::st_geometry_type(habitat))
  bad <- !gt %in% c("POLYGON", "MULTIPOLYGON")
  if (any(bad)) {
    message("  dropping ", sum(bad), " non-polygon feature(s): ",
            paste(unique(gt[bad]), collapse = ", "))
    habitat <- habitat[!bad, ]
  }

  habitat$NewHab <- substr(habitat$NewHabStrat, 1, 2)
  found   <- sort(unique(habitat$NewHab))
  missing <- setdiff(hab.classes, found)
  extra   <- setdiff(found, hab.classes)
  if (length(missing))
    warning("Classes requested but absent from the gdb: ",
            paste(missing, collapse = ", "))
  if (length(extra))
    message("  ignoring unrequested class(es): ", paste(extra, collapse = ", "))

  # ---- points ---------------------------------------------------------------
  # The microgrid table already carries X/Y as EPSG:4326 (verified to match the
  # geometry centroid exactly), so the footprint needs no reprojection at all.
  stopifnot(all(c("X", "Y", "Shape_Area") %in% names(microgrid)))
  mg.pts <- terra::vect(data.frame(x = microgrid$X, y = microgrid$Y,
                                   area = microgrid$Shape_Area),
                        geom = c("x", "y"), crs = "EPSG:4326")

  # Habitat centroids are taken in the native projected CRS (Florida GDL
  # Albers) so they are true centroids, then carried to lon/lat.
  message("Computing habitat centroids...")
  ctr  <- sf::st_transform(sf::st_centroid(sf::st_geometry(habitat)), 4326)
  cxy  <- sf::st_coordinates(ctr)
  hab.pts <- terra::vect(data.frame(x = cxy[, 1], y = cxy[, 2],
                                    area  = habitat$Shape_Area,
                                    class = habitat$NewHab),
                         geom = c("x", "y"), crs = "EPSG:4326")

  # ---- optional QC ----------------------------------------------------------
  if (isTRUE(run.checks))
    fn.check_GFISHER_join(habitat, microgrid, dir.maps)

  # ---- mapping footprint ----------------------------------------------------
  message("Rasterizing mapping footprint...")
  microgrid.ras <- terra::rasterize(mg.pts, depth, field = "area", fun = "sum")
  microgrid.ras[is.na(depth)] <- NA
  names(microgrid.ras) <- "mapped_area_m2"
  mapped <- !is.na(microgrid.ras)
  message(sprintf("  %d of %d water cells have side-scan coverage (%.1f%%)",
                  sum(terra::values(mapped, mat = FALSE), na.rm = TRUE),
                  sum(!is.na(terra::values(depth, mat = FALSE))),
                  100 * sum(terra::values(mapped, mat = FALSE), na.rm = TRUE) /
                        sum(!is.na(terra::values(depth, mat = FALSE)))))

  # ---- per-class proportions ------------------------------------------------
  hab.list <- vector("list", length(hab.classes))
  names(hab.list) <- hab.classes
  for (cls in hab.classes) {
    #cls="NH"
    message("Class ", cls, "...")
    pts.i <- hab.pts[hab.pts$class == cls, ]

    if (nrow(pts.i) == 0L) {
      prop.i <- terra::init(terra::rast(depth), NA)
    } else {
      area.i <- terra::rasterize(pts.i, depth, field = "area", fun = "sum",
                                 background = 0)
      # Proportion of the *surveyed* area in each cell that is this class.
      # microgrid.ras is NA outside the footprint, so prop.i is NA there too.
      prop.i <- area.i / microgrid.ras
    }
    prop.i[!mapped] <- NA

    # Digitised habitat can slightly exceed the scanned area it sits in.
    n.over <- sum(terra::values(prop.i, mat = FALSE) > 1, na.rm = TRUE)
    if (n.over > 0) {
      message("  clamping ", n.over, " cell(s) with proportion > 1")
      prop.i <- terra::clamp(prop.i, 0, 1, values = TRUE)
    }

    if (cls %in% extrapolate && fill.method != "none") {
      attr(prop.i, "mapped.area") <- terra::values(microgrid.ras, mat = FALSE)
      prop.i <- fn.fill_habitat_gaps(prop.i, depth = depth, mapped = mapped,
                                     max.depth = max.depth.hab,
                                     method = fill.method,
                                     covariates = covariates,
                                     strata = strata,
                                     strata.breaks = strata.breaks,
                                     idw.idp = idw.idp,
                                     idw.nmax = idw.nmax,
                                     anchor.zero = anchor.zero,
                                     clamp.to.observed = clamp.to.observed)
    }
    plot(prop.i,main=cls)
    # Anything still unfilled is water with no evidence of this habitat.
    prop.i[is.na(prop.i) & !is.na(depth)] <- 0
    # Land is NA in every layer -- artificial classes included. The previous
    # version left land as 0 for AL/AM/AH.
    prop.i[is.na(depth)] <- NA

    names(prop.i) <- cls
    hab.list[[cls]] <- prop.i
    #plot(prop.i,main=cls)
    rm(pts.i); gc(verbose = FALSE)
  }
  hab.stack <- terra::rast(hab.list)
  names(hab.stack) <- hab.classes
  #plot(hab.stack)
  
  # ---- cap the reef total ---------------------------------------------------
  # Classes are extrapolated independently, so their predictions can sum past
  # 1 in a cell. Left alone that would leave no room for seagrass and seabed in
  # the downstream sum-to-1 step, so rescale the offending cells proportionally.
  if (isTRUE(cap.total)) {
    tot <- sum(hab.stack)
    n.over <- sum(terra::values(tot, mat = FALSE) > 1, na.rm = TRUE)
    if (n.over > 0) {
      message(sprintf("Rescaling %d cell(s) where the reef total exceeded 1 (max %.3f).",
                      n.over, max(terra::values(tot, mat = FALSE), na.rm = TRUE)))
      scale.r <- terra::ifel(tot > 1, 1 / tot, 1)
      hab.stack <- hab.stack * scale.r
      names(hab.stack) <- hab.classes
    }
  }

  if (isTRUE(do.plot)) terra::plot(hab.stack, colNA = "gray90")

  # ---- write ----------------------------------------------------------------
  if (isTRUE(write.ascii)) {
    for (cls in hab.classes) {
      f <- file.path(dir.maps, sprintf("GFISHER_%s_prop_%dmin.asc", cls, res.min))
      terra::writeRaster(hab.stack[[cls]], f, overwrite = TRUE, NAflag = -9999,
                         gdal = c("DECIMAL_PRECISION=6"))
      unlink(paste0(f, ".aux.xml"))
    }
    f <- file.path(dir.maps, sprintf("GFISHER_microgrid_%dmin.asc", res.min))
    terra::writeRaster(microgrid.ras, f, overwrite = TRUE, NAflag = -9999,
                       gdal = c("DECIMAL_PRECISION=2"))
    unlink(paste0(f, ".aux.xml"))
    message("Ecospace ascii files written to\n  ", dir.maps)
  }

  invisible(list(habitat = hab.stack, microgrid = microgrid.ras))
}


#' QC the habitat-polygon to microgrid join.
#'
#' Diagnostic only -- nothing in the rasterisation depends on the join. Split
#' out of the main function because the intersect is expensive.
#'
#' @param habitat sf of habitat polygons.
#' @param microgrid sf of mapped microgrids.
#' @param dir.maps Directory for the exception CSV.
#' @return Invisible data.frame of microgrids whose habitat area exceeds their
#'   grid area.
fn.check_GFISHER_join <- function(habitat, microgrid, dir.maps) {

  message("Joining habitat centroids to microgrids (QC)...")
  ctr <- sf::st_centroid(sf::st_geometry(habitat))
  hit <- sf::st_intersects(ctr, sf::st_geometry(microgrid))

  # The old code did unlist(hit), which silently misaligns as soon as any
  # centroid matches zero or more than one microgrid. Take the first match per
  # centroid and keep the vector the same length as habitat.
  n   <- lengths(hit)
  idx <- rep(NA_integer_, length(hit))
  idx[n > 0] <- vapply(hit[n > 0], function(i) i[1], integer(1))

  message(sprintf("  centroids matching 0 microgrids: %d", sum(n == 0)))
  message(sprintf("  centroids matching >1 microgrid: %d", sum(n > 1)))

  hab <- sf::st_drop_geometry(habitat)
  mcg <- sf::st_drop_geometry(microgrid)
  hab$MicroGrid <- mcg$MicroGrid[idx]
  names(hab)[names(hab) == "Shape_Area"] <- "Hab_Area"
  names(mcg)[names(mcg) == "Shape_Area"] <- "Grid_Area"

  message(sprintf("  habitat records with a matching microgrid: %d of %d",
                  sum(!is.na(hab$MicroGrid)), nrow(hab)))
  message(sprintf("  microgrids with no habitat record: %d",
                  sum(!mcg$MicroGrid %in% hab$MicroGrid)))

  hab.sum <- stats::aggregate(Hab_Area ~ MicroGrid, data = hab, FUN = sum)
  chk <- merge(hab.sum, mcg[, c("MicroGrid", "Grid_Area")],
               by = "MicroGrid", all = TRUE)
  chk$habpct <- round(chk$Hab_Area / chk$Grid_Area, 4)
  over <- chk[which(chk$habpct > 1), ]
  message(sprintf("  microgrids with Hab_Area > Grid_Area: %d", nrow(over)))
  if (nrow(over) > 0)
    utils::write.csv(over, file.path(dir.maps, "hab area vs grid area.csv"),
                     row.names = FALSE)

  invisible(over)
}


#' Colour scale for GFISHER habitat proportion maps.
#'
#' One definition, shared by every GFISHER figure so panels are comparable.
#'
#' Exact zero gets its own bin and a neutral colour; the ramp is spent entirely
#' on cells that actually hold reef. Coverage is so zero-inflated that an even
#' ramp renders the whole shelf as one flat colour, which is what made the
#' first version of these maps unreadable. Breaks come from the quantiles of
#' the positive values.
#'
#' The first break straddles zero rather than ending at it: terra's interval
#' test is left-closed, so a bin ending exactly at 0 pushes the zeros into the
#' next bin and the whole map reads as habitat.
#'
#' @param v Numeric values the scale is derived from. Pass the pooled values of
#'   several rasters to put them all on one comparable scale.
#' @param style "quantile" gives zero its own white bin and spreads the ramp
#'   over the quantiles of the positive values. "pretty" is the original
#'   scheme: evenly spaced breaks over the full range with a bias-2 spline
#'   ramp, so zero takes the low end of the ramp rather than a separate colour.
#' @param n.bins Number of quantile bins ("quantile") or target number of
#'   intervals passed to pretty() ("pretty").
#' @return list(breaks, col, labels). `labels` is NULL when terra should label
#'   the legend itself.
fn.gfisher_scale <- function(v, style = c("quantile", "matlab", "pretty"),
                             n.bins = 9) {

  style <- match.arg(style)
  pos <- v[!is.na(v) & v > 0]
  if (length(pos) < 2)
    return(list(breaks = NULL, col = NULL, labels = NULL, type = "interval"))

  if (style == "matlab") {
    # The original look: a smooth matlab.like2 colourbar, bias-2 so the ramp
    # spends most of its range on the small values this data is made of.
    # raster::plot drew this as a continuous bar; terra needs to be told,
    # otherwise 50 pretty() breaks become 50 discrete legend entries.
    return(list(breaks = NULL,
                col    = colorRampPalette(colorRamps::matlab.like2(100),
                                          bias = 2, interpolate = "spline")(100),
                labels = NULL,
                type   = "continuous"))
  }

  if (style == "pretty") {
    brks <- unique(sort(c(0, min(pos),
                          utils::tail(pretty(v[!is.na(v)], n = 50), -1))))
    brks <- brks[is.finite(brks)]
    if (max(brks) < max(pos)) brks <- c(brks, max(pos))
    # Zero lands in the first interval [0, min(pos)) under terra's left-closed
    # test, so it takes the low end of the ramp -- the original look.
    cols <- colorRampPalette(colorRamps::matlab.like2(length(brks)),
                             bias = 2, interpolate = "spline")(length(brks) - 1)
    return(list(breaks = brks, col = cols, labels = NULL))
  }

  qb <- unique(stats::quantile(pos, seq(0, 1, length.out = n.bins), na.rm = TRUE))
  qb <- qb[qb > 1e-12]

  brks <- c(-1, 1e-12, qb)
  cols <- c("white", colorRamps::matlab.like2(length(brks) - 2))
  lo   <- c(0, utils::head(qb, -1))

  list(breaks = brks,
       col    = cols,
       labels = c("0", paste(signif(lo, 2), "-", signif(qb, 2))))
}


#' Compare gap-fill methods side by side.
#'
#' One row per habitat class, one column per method, sharing a colour scale
#' within each row so the methods are directly comparable. All methods are
#' identical inside the side-scan footprint by construction -- the differences
#' are entirely in the cells nobody surveyed, which is what this figure is for.
#'
#' @param runs Named list of habitat SpatRasters, one per method, as returned
#'   in the `habitat` element of fn.make_GFISHER_habitat_maps().
#' @param dir.maps Output directory for the PNG.
#' @param classes Classes to show, one row each.
#' @param microgrid Optional footprint raster; its outline is drawn on each panel.
#' @param class.labels Optional row labels.
#' @return The path written, invisibly.
fn.plot_GFISHER_fill_comparison <- function(runs, dir.maps,
                                            classes = c("NL", "NM", "NH"),
                                            microgrid = NULL,
                                            class.labels = NULL,
                                            style = c("quantile", "pretty")) {

  style <- match.arg(style)

  stopifnot(length(runs) > 0, !is.null(names(runs)))
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)

  res.min <- round(terra::res(runs[[1]])[1] * 60, 0)
  if (is.null(class.labels))
    class.labels <- c(AL = "Artificial, low relief",  AM = "Artificial, medium relief",
                      AH = "Artificial, high relief", NL = "Natural, low relief",
                      NM = "Natural, medium relief",  NH = "Natural, high relief")[classes]

  fp <- NULL
  if (!is.null(microgrid))
    fp <- try(terra::as.polygons(!is.na(microgrid), dissolve = TRUE), silent = TRUE)
  if (inherits(fp, "try-error")) fp <- NULL

  fout <- file.path(dir.maps, sprintf("GFISHER fill method comparison %dmin %s.png",
                                      res.min, style))
  png(fout, width = 4.3 * length(runs), height = 4 * length(classes),
      units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(length(classes), length(runs)), oma = c(0, 2.5, 3, 0))

  at.row <- seq(1, 0, length.out = 2 * length(classes) + 1)[seq(2, by = 2,
                                                                length.out = length(classes))]
  for (ci in seq_along(classes)) {
    cl <- classes[ci]
    # Pool across methods so one scale serves the whole row.
    s <- fn.gfisher_scale(unlist(lapply(runs, function(r)
      terra::values(r[[cl]], mat = FALSE))), style = style)

    for (mi in seq_along(runs)) {
      r <- runs[[mi]][[cl]]
      v <- terra::values(r, mat = FALSE)
      last <- mi == length(runs)
      terra::plot(r, col = s$col, breaks = s$breaks, colNA = "gray85",
                  mar = c(2, 2, 3.5, if (last) 6 else 1),
                  legend = last,
                  plg = if (!last) NULL
                        else if (is.null(s$labels)) list(cex = 0.85)
                        else list(cex = 0.85, legend = s$labels),
                  main = sprintf("%s\n%d cells > 0   mean %.4f   max %.3f",
                                 names(runs)[mi], sum(v > 0, na.rm = TRUE),
                                 mean(v, na.rm = TRUE), max(v, na.rm = TRUE)))
      if (!is.null(fp)) terra::lines(fp, col = "gray30", lwd = 0.6)
      maps::map("state", add = TRUE, col = "gray40")
    }
    mtext(class.labels[ci], side = 2, outer = TRUE, line = 0.5,
          font = 2, cex = 0.9, at = at.row[ci])
  }
  mtext(sprintf("GFISHER gap-fill comparison (%d min; outline = side-scan footprint)",
                res.min), side = 3, outer = TRUE, line = 0.6, font = 2, cex = 1.05)

  invisible(fout)
}


#' Plot the GFISHER habitat proportion maps and the mapping footprint.
#'
#' @param x Either the list returned by fn.make_GFISHER_habitat_maps(), a
#'   SpatRaster of the habitat classes, or a directory of written ASCII grids.
#' @param dir.maps Output directory for the PNGs.
#' @param hab.labels Optional panel labels, one per class.
#' @param tag Optional suffix for the filenames, e.g. "unfilled".
#' @param col Colour ramp passed straight to terra::plot(). No breaks are
#'   supplied, so terra draws a continuous colourbar over each layer's range.
#'   Note that coverage is heavily zero-inflated, so on a linear ramp most of
#'   the shelf sits at the low end and only the few high cells pick up colour.
#' @return Invisible character vector of the files written.
fn.plot_GFISHER_habitats <- function(x, dir.maps, hab.labels = NULL, tag = "",
                                     col = colorRamps::matlab.like(100)) {

  microgrid <- NULL

  if (is.character(x)) {
    files <- list.files(x, pattern = "\\.asc$", full.names = TRUE)
    is.mg <- grepl("microgrid", basename(files), ignore.case = TRUE)
    hab   <- terra::rast(files[!is.mg])
    # "GFISHER_AH_prop_5min.asc" -> "AH"
    names(hab) <- vapply(strsplit(basename(files[!is.mg]), "_"),
                         `[`, character(1), 2)
    if (any(is.mg)) microgrid <- terra::rast(files[is.mg][1])
  } else if (inherits(x, "SpatRaster")) {
    hab <- x
  } else {
    hab <- x$habitat
    microgrid <- x$microgrid
  }

  if (missing(dir.maps)) stop("dir.maps is required.")
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)

  res.min <- round(terra::res(hab)[1] * 60, 0)
  ord <- intersect(c('AL','AM','AH','NL','NM','NH'), names(hab))
  hab <- hab[[ord]]

  if (is.null(hab.labels))
    hab.labels <- c(AL = "Artificial, low relief",
                    AM = "Artificial, medium relief",
                    AH = "Artificial, high relief",
                    NL = "Natural, low relief",
                    NM = "Natural, medium relief",
                    NH = "Natural, high relief")[ord]

  out <- character(0)

  # Outline of the side-scan footprint, so a cell that is zero because it was
  # surveyed and had no reef is distinguishable from one that is zero because
  # nobody looked.
  fp <- NULL
  if (!is.null(microgrid))
    fp <- try(terra::as.polygons(!is.na(microgrid), dissolve = TRUE), silent = TRUE)
  if (inherits(fp, "try-error")) fp <- NULL

  f1 <- file.path(dir.maps, sprintf("GFISHER habitat prop %dmin%s.png",
                                    res.min,
                                    if (nzchar(tag)) paste0(" ", tag) else ""))
  png(f1, height = 8, width = 8.5, units = "in", res = 300)
  par(mfcol = c(3, 2))
  for (i in seq_len(terra::nlyr(hab))) {
    r.i <- hab[[i]]
    v   <- terra::values(r.i, mat = FALSE)
    pos <- v[!is.na(v) & v > 0]

    # Plain continuous terra::plot -- no breaks, so terra draws its own
    # colourbar over the layer's range.
    terra::plot(r.i, col = col, colNA = "gray85",
                mar = c(2, 2, 3.5, 5.5), plg = list(cex = 0.7),
                main = if (length(pos) > 0)
                         sprintf("%s\n%d cells > 0   max %.3f",
                                 hab.labels[i], length(pos), max(pos))
                       else paste(hab.labels[i], "\n(no coverage)"))
    if (!is.null(fp)) terra::lines(fp, col = "gray30", lwd = 0.6)
    maps::map("state", add = TRUE, col = "gray40")
  }
  dev.off()
  out <- c(out, f1)

  # Mapping footprint, km2 per cell.
  if (!is.null(microgrid)) {
    f2 <- file.path(dir.maps, sprintf("GFISHER mapping footprint %dmin%s.png",
                                      res.min, if (nzchar(tag)) paste0(" ", tag) else ""))
    png(f2, height = 7, width = 7, units = "in", res = 300)
    mg <- microgrid / 1e6
    terra::plot(mg, col = colorRamps::matlab.like2(50), colNA = "gray90",
                mar = c(3, 3, 3, 6), plg = list(cex = 0.7),
                main = "Side-scan mapping footprint\n(sq km per grid cell)")
    maps::map("state", add = TRUE, col = "gray40")
    dev.off()
    out <- c(out, f2)
  }

  message("Figures written to\n  ", dir.maps)
  invisible(out)
}
