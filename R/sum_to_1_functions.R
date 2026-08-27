#' Combine the habitat layers into an Ecospace basemap that sums to 1.
#'
#' Rewrite of "maps/scale habitats to sum1.R" with three changes:
#'   - dbSEABED rock is dissolved into the natural reef classes rather than
#'     kept as a substrate fraction,
#'   - the artificial reef database is merged into the GFISHER artificial
#'     classes,
#'   - the legacy relief multipliers (AH x 4, AM x 2) are dropped, because
#'     fn.make_AR_maps() already weights AR footprint by relief height.
#'
#' Output layers, summing to 1 in every water cell:
#'   AL AM AH        artificial reef, low / medium / high relief
#'   NL NM NH        natural reef,    low / medium / high relief
#'   seagrass
#'   Sand Mud Gravel
#'
#' terra only. Do not add library(raster) or library(sp) here -- see the note at
#' the top of R/GFISHER functions.R.
#'
#' Caveat worth carrying forward: rock exposure averages ~0.114 across the shelf
#' against GFISHER natural reef of ~0.016, so once dissolved roughly 90% of the
#' natural reef in the finished basemap is rock-derived rather than observed by
#' side-scan.


#' Build a depth x latitude stratification on the model grid.
#'
#' @param depth SpatRaster of positive depths (the template grid).
#' @param depth.breaks Depth bin edges, metres.
#' @param lat.breaks Latitude bin edges, decimal degrees.
#' @return list(strata = SpatRaster of integer stratum IDs,
#'   key = data.frame mapping ID to its depth and latitude bin,
#'   depth.bin = SpatRaster of depth-bin IDs, used for the fallback).
fn.habitat_strata <- function(depth,
                              depth.breaks = c(0, 20, 40, 60, 100, 200, Inf),
                              lat.breaks   = seq(25, 31, 1)) {

  dv <- terra::values(depth, mat = FALSE)
  lat <- terra::yFromCell(depth, seq_len(terra::ncell(depth)))

  dbin <- cut(dv,  depth.breaks, include.lowest = TRUE)
  lbin <- cut(lat, lat.breaks,   include.lowest = TRUE)

  key <- expand.grid(dbin = levels(dbin), lbin = levels(lbin),
                     stringsAsFactors = FALSE)
  key$id <- seq_len(nrow(key))

  id <- key$id[match(paste(dbin, lbin), paste(key$dbin, key$lbin))]
  id[is.na(dv)] <- NA

  s <- terra::rast(depth)
  terra::values(s) <- id
  names(s) <- "stratum"

  db <- terra::rast(depth)
  terra::values(db) <- as.integer(dbin)
  names(db) <- "depth_bin"

  list(strata = s, key = key, depth.bin = db)
}


#' Estimate the NL:NM:NH composition of natural reef per stratum.
#'
#' Shares are taken from **mapped cells only**. Unmapped cells hold extrapolated
#' values, and using them would feed the gap-fill model's own output back into
#' the composition it is being used to build.
#'
#' Coverage is uneven -- at 5 min, of 36 depth x latitude strata roughly half
#' have 30 or more cells with reef and a handful have none -- so a stratum with
#' too little data falls back to its depth bin, and then to the whole shelf.
#'
#' @param gfisher SpatRaster with layers NL, NM, NH.
#' @param mapped Logical SpatRaster, TRUE where side-scan mapping exists.
#' @param strata Output of fn.habitat_strata().
#' @param classes Natural classes, in order.
#' @param min.n Minimum cells with reef before a stratum's own shares are used.
#' @param verbose Report the fallback tally.
#' @return list(shares = 3-layer SpatRaster summing to 1,
#'   table = data.frame of the per-stratum shares and source used).
fn.rock_relief_split <- function(gfisher, mapped, strata,
                                 classes = c("NL", "NM", "NH"),
                                 min.n = 20, verbose = TRUE) {

  stopifnot(all(classes %in% names(gfisher)))

  d <- data.frame(
    cell    = seq_len(terra::ncell(gfisher)),
    stratum = terra::values(strata$strata,    mat = FALSE),
    dbin    = terra::values(strata$depth.bin, mat = FALSE),
    mapped  = as.logical(terra::values(mapped, mat = FALSE)))
  for (cl in classes) d[[cl]] <- terra::values(gfisher[[cl]], mat = FALSE)
  d$tot <- rowSums(d[, classes], na.rm = TRUE)

  obs <- d[d$mapped & !is.na(d$tot) & d$tot > 0, , drop = FALSE]
  if (nrow(obs) == 0) stop("No mapped cells with natural reef; cannot split rock.")

  share_of <- function(z) {
    s <- vapply(classes, function(cl) sum(z[[cl]], na.rm = TRUE), numeric(1))
    if (sum(s) <= 0) return(NULL)
    s / sum(s)
  }

  global <- share_of(obs)
  by.dbin <- lapply(split(obs, obs$dbin), share_of)
  by.strat <- lapply(split(obs, obs$stratum), share_of)
  n.strat <- table(obs$stratum)
  n.dbin  <- table(obs$dbin)

  # Resolve one row per stratum present on the grid, recording which level of
  # the fallback chain actually supplied the numbers.
  ids <- sort(unique(stats::na.omit(d$stratum)))
  tab <- do.call(rbind, lapply(ids, function(id) {
    db <- unique(stats::na.omit(d$dbin[d$stratum == id]))[1]
    n  <- if (as.character(id) %in% names(n.strat)) as.integer(n.strat[as.character(id)]) else 0L
    sh <- by.strat[[as.character(id)]]
    src <- "stratum"
    if (is.null(sh) || n < min.n) {
      sh <- by.dbin[[as.character(db)]]
      src <- "depth bin"
      if (is.null(sh) || (as.character(db) %in% names(n.dbin) &&
                          n.dbin[as.character(db)] < min.n)) {
        sh <- global
        src <- "global"
      }
    }
    data.frame(stratum = id, depth_bin = db, n_cells = n, source = src,
               setNames(as.list(round(sh, 4)), classes),
               stringsAsFactors = FALSE)
  }))

  if (verbose) {
    message("  relief composition source: ",
            paste(sprintf("%s=%d", names(table(tab$source)), table(tab$source)),
                  collapse = "  "))
    message(sprintf("  global shares: %s",
                    paste(sprintf("%s %.3f", classes, global), collapse = "  ")))
  }

  # Paint the shares back onto the grid.
  sv <- terra::values(strata$strata, mat = FALSE)
  lyrs <- lapply(classes, function(cl) {
    r <- terra::rast(strata$strata)
    terra::values(r) <- tab[[cl]][match(sv, tab$stratum)]
    r
  })
  shares <- terra::rast(lyrs)
  names(shares) <- classes

  list(shares = shares, table = tab)
}


#' Fill NA gaps in a raster with a focal mean, then zero.
#'
#' dbSEABED has scattered holes; the legacy script smoothed with a 3x3 focal
#' mean before use. Any cell still NA inside the domain becomes 0.
#'
#' @param r SpatRaster.
#' @param depth Template; cells that are NA here stay NA.
#' @param w Focal window size.
#' @return SpatRaster with no NA inside the domain.
fn.fill_na_focal <- function(r, depth, w = 3) {
  out <- r
  if (any(is.na(terra::values(out, mat = FALSE)) &
          !is.na(terra::values(depth, mat = FALSE)))) {
    sm <- terra::focal(out, w = w, fun = "mean", na.rm = TRUE, na.policy = "only")
    out <- terra::cover(out, sm)
  }
  out[is.na(out) & !is.na(depth)] <- 0
  out[is.na(depth)] <- NA
  out
}


#' Combine every habitat source into layers that sum to 1.
#'
#' @param depth SpatRaster of positive depths (the template grid).
#' @param gfisher SpatRaster with AL, AM, AH, NL, NM, NH (proportions).
#' @param ar SpatRaster from fn.make_AR_maps(); layers Low, Medium, High are
#'   used, "all" is ignored.
#' @param seabed SpatRaster from fn.rasterize_dbseabed(); layer names are matched
#'   case-insensitively for RCK, SND, MUD, GVL.
#' @param seagrass Single-layer SpatRaster of seagrass cover proportion.
#' @param microgrid SpatRaster of scanned area, used to identify mapped cells.
#' @param dir.maps Output directory.
#' @param depth.breaks,lat.breaks Stratification for the rock relief split.
#' @param min.n Minimum cells per stratum before falling back.
#' @param dissolve.rock Dissolve rock into the natural classes. FALSE keeps RCK
#'   as an eleventh layer in the sum instead.
#' @param write.ascii Write Ecospace ASCII grids.
#' @param tol Tolerance for the sum-to-1 assertion.
#' @param verbose Print progress and QC.
#' @return list(habitat, rock_raw, qc, shares).
fn.combine_habitats_sum1 <- function(depth, gfisher, ar, seabed, seagrass,
                                     microgrid,
                                     dir.maps,
                                     depth.breaks  = c(0, 20, 40, 60, 100, 200, Inf),
                                     lat.breaks    = seq(25, 31, 1),
                                     min.n         = 20,
                                     dissolve.rock = TRUE,
                                     write.ascii   = TRUE,
                                     tol           = 1e-6,
                                     verbose       = TRUE) {

  if (inherits(depth, "Raster")) depth <- terra::rast(depth)
  if (is.na(terra::crs(depth)) || !nzchar(terra::crs(depth)))
    terra::crs(depth) <- "EPSG:4326"
  res.min <- round(terra::res(depth)[1] * 60, 0)
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)

  water <- !is.na(terra::values(depth, mat = FALSE))
  zero  <- terra::ifel(is.na(depth), NA, 0)

  # ---- substrate ------------------------------------------------------------
  pick <- function(x, pat) {
    i <- grep(pat, names(x), ignore.case = TRUE)
    if (length(i) == 0) stop("No layer matching '", pat, "' in: ",
                             paste(names(x), collapse = ", "))
    x[[i[1]]]
  }
  rock   <- fn.fill_na_focal(pick(seabed, "RCK"), depth)
  sand   <- fn.fill_na_focal(pick(seabed, "SND"), depth)
  mud    <- fn.fill_na_focal(pick(seabed, "MUD"), depth)
  gravel <- fn.fill_na_focal(pick(seabed, "GVL"), depth)
  rock.raw <- rock
  names(rock.raw) <- "RCK_raw"

  # ---- 1. artificial: GFISHER + the AR database, matched by name ------------
  if (verbose) message("Combining artificial reef classes...")
  ar.map <- c(AL = "Low", AM = "Medium", AH = "High")
  A <- list()
  for (cl in names(ar.map)) {
    a <- gfisher[[cl]]
    j <- which(tolower(names(ar)) == tolower(ar.map[[cl]]))
    if (length(j) == 1) {
      a <- a + terra::ifel(is.na(ar[[j]]), 0, ar[[j]])
    } else if (verbose) {
      message("    no AR layer '", ar.map[[cl]], "'; using GFISHER alone for ", cl)
    }
    a[is.na(depth)] <- NA
    A[[cl]] <- a
  }

  # ---- 2. natural: GFISHER + rock split by stratum --------------------------
  strata <- fn.habitat_strata(depth, depth.breaks, lat.breaks)
  mapped <- !is.na(microgrid)
  split  <- fn.rock_relief_split(gfisher, mapped, strata, min.n = min.n,
                                 verbose = verbose)

  N <- list()
  for (cl in c("NL", "NM", "NH")) {
    n <- gfisher[[cl]]
    if (isTRUE(dissolve.rock)) {
      add <- rock * split$shares[[cl]]
      add[is.na(add)] <- 0
      n <- n + add
    }
    n[is.na(depth)] <- NA
    N[[cl]] <- n
  }

  # ---- 3. cap the reef total ------------------------------------------------
  reef <- terra::rast(c(A, N))
  names(reef) <- c(names(A), names(N))
  if (!isTRUE(dissolve.rock)) reef <- c(reef, setNames(rock, "RCK"))

  tot.reef <- sum(reef)
  n.over <- sum(terra::values(tot.reef, mat = FALSE) > 1, na.rm = TRUE)
  if (n.over > 0) {
    if (verbose)
      message(sprintf("  rescaling %d cell(s) where reef + rock exceeded 1 (max %.3f)",
                      n.over, max(terra::values(tot.reef, mat = FALSE), na.rm = TRUE)))
    sc <- terra::ifel(tot.reef > 1, 1 / tot.reef, 1)
    nm <- names(reef)
    reef <- reef * sc
    names(reef) <- nm
    tot.reef <- sum(reef)
  }

  # ---- 4. seagrass takes what it can from the remainder ---------------------
  remain <- 1 - tot.reef
  remain[remain < 0] <- 0
  sav <- terra::ifel(is.na(seagrass), 0, seagrass)
  sav <- min(sav, remain)
  sav[is.na(depth)] <- NA
  names(sav) <- "seagrass"
  remain <- remain - sav

  # ---- 5. sediment fills the rest, renormalised among itself ----------------
  sed <- terra::rast(list(Sand = sand, Mud = mud, Gravel = gravel))
  names(sed) <- c("Sand", "Mud", "Gravel")
  sed.sum <- sum(sed)

  # Where dbSEABED has nothing to say, put the remainder in Sand so the cell
  # still closes to 1 rather than leaving a hole.
  no.sed <- sed.sum <= 0
  n.nosed <- sum(terra::values(no.sed, mat = FALSE) & water, na.rm = TRUE)
  if (n.nosed > 0 && verbose)
    message(sprintf("  %d cell(s) have no dbSEABED sediment; remainder assigned to Sand",
                    n.nosed))
  sed[[1]] <- terra::ifel(no.sed, 1, sed[[1]])
  sed.sum  <- sum(sed)
  sed <- (sed / sed.sum) * remain
  for (i in 1:3) sed[[i]][is.na(depth)] <- NA
  names(sed) <- c("Sand", "Mud", "Gravel")

  # ---- assemble + QC --------------------------------------------------------
  habitat <- c(reef, sav, sed)
  total   <- sum(habitat)
  tv <- terra::values(total, mat = FALSE)
  bad <- sum(abs(tv[water] - 1) > tol, na.rm = TRUE)
  if (verbose)
    message(sprintf("  sum check: %d of %d water cells off 1 by more than %g",
                    bad, sum(water), tol))
  if (bad > 0)
    warning(bad, " water cell(s) do not sum to 1 within ", tol, ".")

  qc <- do.call(rbind, lapply(names(habitat), function(n) {
    v <- terra::values(habitat[[n]], mat = FALSE)
    data.frame(layer = n,
               mean_prop = round(mean(v, na.rm = TRUE), 5),
               max_prop  = round(max(v, na.rm = TRUE), 4),
               cells_gt_0 = sum(v > 0, na.rm = TRUE),
               cells_gt_10pct = sum(v > 0.1, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  if (verbose) print(qc, row.names = FALSE)

  # ---- write ----------------------------------------------------------------
  if (isTRUE(write.ascii)) {
    for (n in names(habitat)) {
      f <- file.path(dir.maps, sprintf("habitat_%s_%dmin.asc", n, res.min))
      terra::writeRaster(habitat[[n]], f, overwrite = TRUE, NAflag = -9999,
                         gdal = c("DECIMAL_PRECISION=6"))
      unlink(paste0(f, ".aux.xml"))
    }
    f <- file.path(dir.maps, sprintf("RCK_raw_%dmin.asc", res.min))
    terra::writeRaster(rock.raw, f, overwrite = TRUE, NAflag = -9999,
                       gdal = c("DECIMAL_PRECISION=6"))
    unlink(paste0(f, ".aux.xml"))

    utils::write.csv(qc, file.path(dir.maps,
                                   sprintf("habitat_basemap_QC_%dmin.csv", res.min)),
                     row.names = FALSE)
    utils::write.csv(split$table,
                     file.path(dir.maps,
                               sprintf("rock_relief_shares_%dmin.csv", res.min)),
                     row.names = FALSE)
    if (verbose) message("Ecospace ascii files written to\n  ", dir.maps)
  }

  invisible(list(habitat = habitat, rock_raw = rock.raw, qc = qc,
                 shares = split$table, total = total))
}


#' Plot the finished habitat basemap.
#'
#' @param x list from fn.combine_habitats_sum1(), or a SpatRaster of the layers.
#' @param dir.maps Output directory for the PNG.
#' @param col Colour ramp passed to terra::plot().
#' @return The path written, invisibly.
fn.plot_habitat_basemap <- function(x, dir.maps,
                                    col = colorRampPalette(
                                      colorRamps::matlab.like2(100),
                                      bias = 2, interpolate = "spline")(100)) {

  hab <- if (inherits(x, "SpatRaster")) x else x$habitat
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)
  res.min <- round(terra::res(hab)[1] * 60, 0)

  stk <- c(hab, setNames(sum(hab), "Sum"))
  fout <- file.path(dir.maps, sprintf("habitat basemap %dmin.png", res.min))

  png(fout, height = 10, width = 9, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(4, 3))
  for (i in seq_len(terra::nlyr(stk))) {
    v <- terra::values(stk[[i]], mat = FALSE)
    terra::plot(stk[[i]], col = col, colNA = "lightgray",
                mar = c(2, 2, 3, 5), plg = list(cex = 0.6),
                main = sprintf("%s\nmean %.4f   max %.3f", names(stk)[i],
                               mean(v, na.rm = TRUE), max(v, na.rm = TRUE)))
    maps::map("state", add = TRUE, fill = TRUE, col = "lightgray")
  }

  message("Figure written to\n  ", fout)
  invisible(fout)
}
