#' Artificial reef basemaps for WFS Ecospace.
#'
#' Rewrite of "maps/Artificial reefs/make AR area map.R" as callable functions.
#'
#' Two sources are combined:
#'   dataS2_artificial_reef_structures_REDACTED.csv  structure records with
#'     coordinates and footprint area (m2), many states; FL_GOM is used here.
#'   reeflocations.csv                               the FWC deployment table,
#'     which carries relief height (m) but is keyed only by free-text
#'     Description, so the two are joined by fuzzy string matching.
#'
#' Pipeline:
#'   1. Fill missing relief within the FWC table by fuzzy-matching descriptions
#'      of records that lack relief against those that have it.
#'   2. Join AR structures to FWC records by description (exact, then
#'      quote-stripped, then fuzzy) to attach a relief height.
#'   3. Split relief into Low / Medium / High by 1-D k-means.
#'   4. Rasterize the per-cell sum onto the depth grid, one layer per class
#'      plus "all", and divide by cell area.
#'
#' terra + base R only. Do not add library(raster) or library(sp) here -- see
#' the note at the top of R/GFISHER functions.R.


#' Best fuzzy match of one string against a vector of candidates.
#'
#' Two-stage: agrep() to shortlist, then adist() to pick the closest. Kept
#' verbatim in behaviour from the original script.
#'
#' @param query Single string.
#' @param candidates Character vector to match against.
#' @param max.distance Passed to agrep(); larger is more permissive.
#' @param ignore.case Passed to agrep().
#' @param normalize Strip punctuation and smart quotes, collapse whitespace.
#' @return list(index, match, distance, n_candidates).
fn.best_agrep_match <- function(query, candidates, max.distance = 0.3,
                                ignore.case = TRUE, normalize = TRUE) {

  stopifnot(is.character(query), length(query) == 1, is.character(candidates))

  norm <- function(x) {
    x <- tolower(x)
    x <- gsub("[‘’]", "'", x)     # curly single quotes
    x <- gsub("[“”]", "\"", x)    # curly double quotes
    x <- gsub("[^[:alnum:] ]", " ", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }
  q  <- if (normalize) norm(query)      else query
  cs <- if (normalize) norm(candidates) else candidates

  hits <- agrep(q, cs, max.distance = max.distance, ignore.case = ignore.case)
  if (length(hits) == 0)
    return(list(index = NA_integer_, match = NA_character_,
                distance = NA_real_, n_candidates = 0))

  d <- adist(q, cs[hits], partial = TRUE, ignore.case = FALSE)
  best <- which.min(d)

  list(index        = hits[best],
       match        = candidates[hits[best]],
       distance     = as.numeric(d[best]),
       n_candidates = length(hits))
}


#' Vectorised wrapper around fn.best_agrep_match().
#'
#' @inheritParams fn.best_agrep_match
#' @param queries Character vector of strings to match.
#' @return data.frame with one row per query.
fn.best_agrep_join <- function(queries, candidates, max.distance = 0.3,
                               ignore.case = TRUE, normalize = TRUE) {
  do.call(rbind, lapply(queries, function(q) {
    r <- fn.best_agrep_match(q, candidates, max.distance, ignore.case, normalize)
    data.frame(query = q, match = r$match, index = r$index,
               distance = r$distance, n_candidates = r$n_candidates,
               stringsAsFactors = FALSE)
  }))
}


#' Fill missing relief heights in the FWC reef table.
#'
#' Records with no relief are matched by description to records that have one,
#' and inherit that value.
#'
#' @param reef data.frame from reeflocations.csv; needs Description and Relief.
#' @param max.distance Fuzzy tolerance. The original used 0.75, which is very
#'   permissive -- inspect the returned n_filled before trusting it.
#' @param verbose Report how many were filled.
#' @return `reef` with Relief filled where a match was found. Rows are
#'   reordered (relief-present first), as in the original.
fn.fill_reef_relief <- function(reef, max.distance = 0.75, verbose = TRUE) {

  strip <- function(x) tolower(gsub("[^[:alnum:]]", "", x))

  withrel <- reef[which(reef$Relief > 0), ]
  norel   <- reef[which(!(reef$Relief > 0)), ]
  if (nrow(norel) == 0) return(reef)

  query      <- unique(strip(norel$Description))
  candidates <- unique(strip(withrel$Description))

  best <- fn.best_agrep_join(queries = query, candidates = candidates,
                             max.distance = max.distance)

  # Map the stripped forms back to the original text so the join can be made
  # on Description as it actually appears in each table.
  best$match.og <- withrel$Description[match(best$match, strip(withrel$Description))]
  best$query.og <- norel$Description[match(best$query,  strip(norel$Description))]
  best$relief   <- withrel$Relief[match(best$match, strip(withrel$Description))]

  norel$Relief <- best$relief[match(norel$Description, best$query.og)]

  if (verbose)
    message(sprintf("  relief filled for %d of %d FWC records lacking it",
                    sum(!is.na(norel$Relief)), nrow(norel)))

  rbind(withrel, norel)
}


#' Attach relief height to AR structure records.
#'
#' Three passes of increasing permissiveness, matching the original: exact
#' description, then with double quotes stripped, then fuzzy.
#'
#' @param ar data.frame of AR structures; needs `description`.
#' @param reef data.frame from fn.fill_reef_relief(); needs Description, Relief.
#' @param max.distance Fuzzy tolerance for the third pass.
#' @param verbose Report match counts per pass.
#' @return `ar` with a `relief` column added.
fn.match_ar_relief <- function(ar, reef, max.distance = 0.5, verbose = TRUE) {

  strip <- function(x) tolower(gsub("[^[:alnum:]]", "", x))

  idx <- match(tolower(ar$description), tolower(reef$Description))
  n1  <- sum(!is.na(idx))

  miss <- is.na(idx)
  idx[miss] <- match(tolower(gsub('\\"', '', ar$description[miss])),
                     tolower(reef$Description))
  n2 <- sum(!is.na(idx)) - n1

  miss <- is.na(idx)
  if (any(miss)) {
    idx[miss] <- vapply(
      lapply(ar$description[miss], function(x)
        agrep(strip(x), strip(reef$Description), max.distance = max.distance)),
      function(h) if (length(h)) h[1] else NA_integer_, integer(1))
  }
  n3 <- sum(!is.na(idx)) - n1 - n2

  ar$relief <- reef$Relief[idx]

  if (verbose) {
    message(sprintf("  relief matched: %d exact, %d quote-stripped, %d fuzzy; %d unmatched",
                    n1, n2, n3, sum(is.na(idx))))
    message(sprintf("  records with a relief value: %d of %d",
                    sum(!is.na(ar$relief)), nrow(ar)))
  }
  ar
}


#' Split relief heights into Low / Medium / High by 1-D k-means.
#'
#' Cut points are the midpoints between the sorted cluster centres, so the
#' classification is reported as explicit breaks rather than cluster labels.
#'
#' @param relief Numeric vector of relief heights (m).
#' @param k Number of classes.
#' @param labels Class labels, low to high.
#' @param seed Passed to set.seed() -- k-means starts are random.
#' @param verbose Print the centres, breaks and class counts.
#' @return list(class = factor, breaks = numeric, centers = numeric).
fn.classify_ar_relief <- function(relief, k = 3,
                                  labels = c("Low", "Medium", "High"),
                                  seed = 42, verbose = TRUE) {

  stopifnot(length(labels) == k)
  ok <- !is.na(relief)
  if (sum(ok) < k) stop("Not enough non-missing relief values to classify.")

  set.seed(seed)
  km      <- stats::kmeans(relief[ok], centers = k)
  centers <- sort(km$centers)
  cuts    <- c(-Inf, (utils::head(centers, -1) + centers[-1]) / 2, Inf)

  cls <- cut(relief, breaks = cuts, labels = labels, include.lowest = TRUE)

  if (verbose) {
    message("  k-means centres: ", paste(round(centers, 2), collapse = ", "))
    message("  class breaks   : ", paste(round(cuts[c(-1, -length(cuts))], 2),
                                         collapse = ", "))
    message(paste(utils::capture.output(print(table(cls, useNA = "ifany"))),
                  collapse = "\n"))
  }
  list(class = cls, breaks = cuts, centers = as.numeric(centers))
}


#' Build artificial reef coverage maps on the Ecospace grid.
#'
#' Entry point called from make_WFS_basemaps.R.
#'
#' @param depth SpatRaster of positive depths (the template grid).
#' @param file.ar Path to dataS2_artificial_reef_structures_REDACTED.csv.
#' @param file.reef Path to reeflocations.csv.
#' @param dir.maps Output directory.
#' @param state Value of the `state` column to keep.
#' @param weight.by.relief Multiply each structure's footprint by its relief
#'   height before summing, as the original did. See the note below -- this
#'   makes the output an index, not a proportion.
#' @param classes Relief classes to write, plus "all" for the ungrouped total.
#' @param fill.relief Fuzzy-fill missing relief in the FWC table first.
#' @param match.distance,fill.distance Fuzzy tolerances for the two join steps.
#' @param seed Passed to the k-means classification.
#' @param write.ascii Write Ecospace ASCII grids.
#' @param verbose Print progress and join diagnostics.
#' @return SpatRaster with one layer per entry of `classes`.
#'
#' @section Units:
#' With `weight.by.relief = TRUE` (the original behaviour) each cell holds
#' `sum(area_m2 * relief_m) / 1e6 / cell_area_km2`. That is not dimensionless
#' and not bounded by 1 -- median FWC relief is 8 m, so it inflates the true
#' footprint fraction by roughly an order of magnitude. The original file names
#' call it "prop_area", but it is a relief-weighted index. Set
#' `weight.by.relief = FALSE` for the genuine fraction of cell area covered.
fn.make_AR_maps <- function(depth, file.ar, file.reef, dir.maps,
                            state            = "FL_GOM",
                            weight.by.relief = TRUE,
                            classes          = c("all", "Low", "Medium", "High"),
                            fill.relief      = TRUE,
                            match.distance   = 0.5,
                            fill.distance    = 0.75,
                            seed             = 42,
                            write.ascii      = TRUE,
                            verbose          = TRUE) {

  if (inherits(depth, "Raster")) depth <- terra::rast(depth)
  stopifnot(inherits(depth, "SpatRaster"))
  if (is.na(terra::crs(depth)) || !nzchar(terra::crs(depth)))
    terra::crs(depth) <- "EPSG:4326"

  res.min <- round(terra::res(depth)[1] * 60, 0)
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)

  # ---- read ---------------------------------------------------------------
  if (verbose) message("Reading artificial reef data...")
  ar <- utils::read.csv(file.ar, row.names = 1, stringsAsFactors = FALSE,
                        allowEscapes = TRUE)
  reef <- utils::read.csv(file.reef, row.names = 1, stringsAsFactors = FALSE,
                          fileEncoding = "windows-1252")

  ar <- ar[ar$state == state & !is.na(ar$lat_dd) & !is.na(ar$long_dd), ]
  if (nrow(ar) == 0) stop("No records left after filtering to state = '", state, "'.")
  if (verbose) message(sprintf("  %d %s structures with coordinates", nrow(ar), state))

  # ---- relief -------------------------------------------------------------
  if (isTRUE(fill.relief)) {
    if (verbose) message("Filling missing relief in the FWC table...")
    reef <- fn.fill_reef_relief(reef, max.distance = fill.distance,
                                verbose = verbose)
  }
  if (verbose) message("Matching structures to FWC relief...")
  ar <- fn.match_ar_relief(ar, reef, max.distance = match.distance,
                           verbose = verbose)

  if (verbose) message("Classifying relief...")
  cl <- fn.classify_ar_relief(ar$relief, seed = seed, verbose = verbose)
  ar$relief_class <- cl$class

  # ---- the quantity being summed ------------------------------------------
  ar$value <- if (isTRUE(weight.by.relief)) ar$area_m2 * ar$relief else ar$area_m2

  n.drop <- sum(is.na(ar$value))
  if (n.drop > 0 && verbose)
    message(sprintf("  %d structure(s) have no %s and contribute nothing",
                    n.drop, if (weight.by.relief) "area or relief" else "area"))

  # ---- rasterize ----------------------------------------------------------
  if (verbose) message("Rasterizing...")
  cell.km2 <- terra::cellSize(depth, unit = "km")

  lyrs <- vector("list", length(classes))
  names(lyrs) <- classes
  for (cls.i in classes) {

    pts.i <- if (cls.i == "all") ar else ar[which(ar$relief_class == cls.i), ]
    pts.i <- pts.i[!is.na(pts.i$value), ]

    if (nrow(pts.i) == 0) {
      r.i <- terra::init(terra::rast(depth), 0)
    } else {
      v <- terra::vect(data.frame(x = pts.i$long_dd, y = pts.i$lat_dd,
                                  value = pts.i$value),
                       geom = c("x", "y"), crs = "EPSG:4326")
      r.i <- terra::rasterize(v, depth, field = "value", fun = "sum",
                              background = 0)
      # km2 of structure (or m3 of structure when relief-weighted) per km2 of cell
      r.i <- (r.i / 1e6) / cell.km2
    }
    r.i[is.na(depth)] <- NA
    names(r.i) <- cls.i
    lyrs[[cls.i]] <- r.i
  }
  ar.stack <- terra::rast(lyrs)
  names(ar.stack) <- classes

  if (verbose) {
    for (cls.i in classes) {
      v <- terra::values(ar.stack[[cls.i]], mat = FALSE)
      message(sprintf("  %-6s cells > 0: %4d   max %.5f   sum %.4f",
                      cls.i, sum(v > 0, na.rm = TRUE), max(v, na.rm = TRUE),
                      sum(v, na.rm = TRUE)))
    }
  }

  # ---- write --------------------------------------------------------------
  if (isTRUE(write.ascii)) {
    for (cls.i in classes) {
      f <- file.path(dir.maps, sprintf("AR_prop_area_%s_%dmin.asc", cls.i, res.min))
      terra::writeRaster(ar.stack[[cls.i]], f, overwrite = TRUE, NAflag = -9999,
                         gdal = c("DECIMAL_PRECISION=8"))
      unlink(paste0(f, ".aux.xml"))
    }
    if (verbose) message("Ecospace ascii files written to\n  ", dir.maps)
  }

  invisible(ar.stack)
}


#' Plot the artificial reef maps.
#'
#' Colour scale carried over from "make AR area map.R": matlab.like2 run
#' through a bias-2 spline ramp, which spends most of the colour range on small
#' values -- necessary here because coverage is heavily zero-inflated.
#'
#' The original passed `breaks = pretty(vals, n = 50)` with
#' `lab.breaks = round(brks, 3)` to raster::plot, which drew a labelled
#' colourbar. terra turns explicit breaks into one discrete legend entry per
#' break, which at n = 50 is unreadable, so the same palette is used against a
#' continuous bar instead. The value-to-colour mapping is effectively identical,
#' since the original breaks were evenly spaced anyway.
#'
#' @param x SpatRaster from fn.make_AR_maps(), or a directory of written ASCII.
#' @param dir.maps Output directory for the PNG.
#' @param tag Optional filename suffix.
#' @param col Colour ramp passed to terra::plot().
#' @return The path written, invisibly.
fn.plot_AR_maps <- function(x, dir.maps, tag = "",
                            col = colorRampPalette(colorRamps::matlab.like2(100),
                                                   bias = 2,
                                                   interpolate = "spline")(100)) {

  if (is.character(x)) {
    files <- list.files(x, pattern = "^AR_prop_area_.*\\.asc$", full.names = TRUE)
    x <- terra::rast(files)
    names(x) <- sub("^AR_prop_area_(.*)_[0-9]+min$", "\\1",
                    tools::file_path_sans_ext(basename(files)))
  }
  if (missing(dir.maps)) stop("dir.maps is required.")
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)

  res.min <- round(terra::res(x)[1] * 60, 0)
  ord <- intersect(c("all", "Low", "Medium", "High"), names(x))
  x   <- x[[ord]]

  fout <- file.path(dir.maps, sprintf("AR prop area %dmin%s.png", res.min,
                                      if (nzchar(tag)) paste0(" ", tag) else ""))
  png(fout, height = 7, width = 7, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  par(mfcol = c(2, 2))

  for (i in seq_len(terra::nlyr(x))) {
    r.i <- x[[i]]
    v   <- terra::values(r.i, mat = FALSE)
    pos <- v[!is.na(v) & v > 0]
    terra::plot(r.i, col = col, colNA = "lightgray",
                mar = c(4, 3, 3, 7), plg = list(cex = 0.7),
                main = if (length(pos) > 0)
                         sprintf("%s\n%d cells > 0   max %.4f",
                                 names(x)[i], length(pos), max(pos))
                       else paste(names(x)[i], "\n(no coverage)"))
    maps::map("state", add = TRUE, fill = TRUE, col = "lightgray")
  }

  message("Figure written to\n  ", fout)
  invisible(fout)
}
