#' Survey region basemaps for WFS Ecospace.
#'
#' Port of "maps/regions/R/*" (itself an R conversion of 27 Python scripts).
#' Builds one categorical grid coding which survey region each model cell
#' belongs to:
#'
#'   -9999  land
#'       0  water, unsampled
#'    1- 9  age-0 survey regions (bays and estuaries, digitized polygons)
#'   10-16  GFISHER survey coverage: 10 FWRI, 11 PASC, 12 PC,
#'          13 F+P, 14 F+PC, 15 P+PC, 16 F+P+PC
#'
#' Where the two overlap, age-0 wins.
#'
#' Not ported: the digitizing step (digitize_georeference_age0.R), which
#' extracts the age-0 polygons from a georeferenced PNG. The legacy README is
#' explicit that R's PNG decoder renders the anti-aliased boundaries about a
#' pixel thinner than Python's, which can leave a gap and resolve only 8 of the
#' 9 regions. `age0_survey_regions.shp` from the Python pipeline is the
#' authoritative input and is treated as source data here.
#'
#' terra + sf only. Note that `sf_use_s2(FALSE)` is set inside the functions
#' that do point-in-polygon: the legacy matched shapely's planar predicate, and
#' spherical geometry gives different edge cells.


#' Gaussian kernel density at the sample points themselves.
#'
#' Reproduces scipy.stats.gaussian_kde with Scott's factor, which is what the
#' original Python used to rank points before trimming to the densest fraction.
#' O(n^2) and deliberately so -- it must match the reference implementation.
#'
#' @param xy Two-column matrix of coordinates.
#' @return Numeric vector of densities, one per row of `xy`.
fn.gaussian_kde_at_points <- function(xy) {
  n <- nrow(xy); d <- 2
  H  <- (n^(-1 / (d + 4)))^2 * stats::cov(xy)   # Scott factor^2 * covariance
  Hi <- solve(H)
  norm <- 1 / (2 * pi * sqrt(det(H)))
  dens <- numeric(n)
  for (k in seq_len(n)) {
    dx <- sweep(xy, 2, xy[k, ], "-")
    dens[k] <- mean(norm * exp(-0.5 * rowSums((dx %*% Hi) * dx)))
  }
  dens
}


#' Cell centres of the depth grid as an sf point layer.
#'
#' terra cell order (row-major from the top-left) is used throughout, so every
#' vector built here lines up with terra::values().
#'
#' @param depth SpatRaster template.
#' @return sf POINT layer, one row per cell.
fn.cell_centres_sf <- function(depth) {
  xy <- terra::xyFromCell(depth, seq_len(terra::ncell(depth)))
  sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]),
               coords = c("x", "y"), crs = 4326)
}


#' Rasterize the age-0 survey region polygons onto the model grid.
#'
#' A cell takes a region ID if its centre falls inside that polygon; other water
#' cells are 0 and land is NA.
#'
#' @param depth SpatRaster template.
#' @param file.shp Path to age0_survey_regions.shp.
#' @param dir.maps Output directory. NULL skips writing.
#' @param region.names Fallback names, used only if the .dbf is missing.
#' @param write.per.region Also write one grid per region.
#' @param verbose Report the cell tally.
#' @return SpatRaster of region IDs.
fn.make_age0_region_grids <- function(depth, file.shp, dir.maps = NULL,
                                      region.names = c("Saint Andrew Bay",
                                        "Saint Joe Bay", "Turkey Point",
                                        "Mid Big Bend", "Cedar Key", "Tampa Bay",
                                        "Sarasota Bay", "Charlotte Harbor",
                                        "Marco Island"),
                                      write.per.region = TRUE,
                                      verbose = TRUE) {

  old.s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)                       # planar, matching shapely
  on.exit(sf::sf_use_s2(old.s2), add = TRUE)

  res.min <- round(terra::res(depth)[1] * 60, 0)

  poly <- sf::st_read(file.shp, quiet = TRUE)
  # Fall back to feature order if the .dbf is not alongside the .shp.
  if (!"ID"   %in% names(poly)) poly$ID   <- seq_len(nrow(poly))
  if (!"NAME" %in% names(poly)) poly$NAME <- region.names[poly$ID]
  poly <- poly[order(poly$ID), ]
  if (verbose) message(sprintf("  %d age-0 region polygons", nrow(poly)))

  pts <- fn.cell_centres_sf(depth)
  hit <- sf::st_intersects(pts, poly)
  regv <- vapply(hit, function(x) if (length(x)) as.integer(poly$ID[x[1]]) else 0L,
                 integer(1))

  r <- terra::rast(depth)
  terra::values(r) <- regv
  r[is.na(depth)] <- NA
  names(r) <- "age0_regions"

  if (verbose) print(table(terra::values(r, mat = FALSE), useNA = "ifany"))

  if (!is.null(dir.maps)) {
    if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)
    fn.write_region_asc(r, file.path(dir.maps,
                                     sprintf("age0_survey_regions_%dmin.asc", res.min)))
    if (isTRUE(write.per.region)) {
      safe <- gsub("[^A-Za-z0-9]+", "", poly$NAME)
      for (k in seq_len(nrow(poly))) {
        id <- poly$ID[k]
        g <- terra::ifel(r == id, id, 0)
        g[is.na(depth)] <- NA
        fn.write_region_asc(g, file.path(dir.maps,
          sprintf("age0_region_%d_%s_%dmin.asc", id, safe[k], res.min)))
      }
    }
  }
  r
}


#' GFISHER survey coverage regions on the model grid.
#'
#' For each survey (FWRI, PASC, PC): rank sample points by 2-D Gaussian kernel
#' density, keep the densest `keep` fraction, take a concave hull, then code
#' each water cell by which hulls contain its centre.
#'
#' @param depth SpatRaster template.
#' @param file.csv GFISHER sample CSV with LAB, LON_DD, LAT_DD.
#' @param dir.maps Output directory. NULL skips writing.
#' @param labs Survey labels, in the order that sets codes 10, 11, 12.
#' @param keep Fraction of points retained by density.
#' @param ratio Passed to sf::st_concave_hull().
#' @param verbose Report points kept per survey and the cell tally.
#' @return SpatRaster of survey codes.
fn.make_gfisher_survey_regions <- function(depth, file.csv, dir.maps = NULL,
                                           labs = c("FWRI", "PASC", "PC"),
                                           keep = 0.95, ratio = 0.10,
                                           verbose = TRUE) {

  old.s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old.s2), add = TRUE)

  res.min <- round(terra::res(depth)[1] * 60, 0)

  d <- utils::read.csv(file.csv)
  need <- c("LAB", "LON_DD", "LAT_DD")
  if (!all(need %in% names(d)))
    stop("GFISHER CSV is missing: ", paste(setdiff(need, names(d)), collapse = ", "))
  d <- d[!is.na(d$LAT_DD) & !is.na(d$LON_DD), need]

  pts <- fn.cell_centres_sf(depth)
  inside <- matrix(FALSE, nrow = nrow(pts), ncol = length(labs),
                   dimnames = list(NULL, labs))

  for (k in labs) {
    xy <- as.matrix(d[d$LAB == k, c("LON_DD", "LAT_DD")])
    if (nrow(xy) == 0) {
      warning("No GFISHER points for survey '", k, "'.")
      next
    }
    dens <- fn.gaussian_kde_at_points(xy)
    sel  <- dens >= stats::quantile(dens, 1 - keep)
    if (verbose) message(sprintf("  %s: kept %d/%d points", k, sum(sel), nrow(xy)))

    mp   <- sf::st_union(sf::st_as_sf(as.data.frame(xy[sel, ]),
                                      coords = c("LON_DD", "LAT_DD"), crs = 4326))
    hull <- sf::st_concave_hull(mp, ratio = ratio, allow_holes = FALSE)
    inside[, k] <- lengths(sf::st_intersects(pts, hull)) > 0
  }

  # Codes are assigned in increasing specificity, so a cell in all three hulls
  # ends on 16 rather than on whichever single survey was written last.
  f <- inside[, 1]; p <- inside[, 2]; c3 <- inside[, 3]
  code <- integer(nrow(pts))
  code[f]            <- 10L
  code[p]            <- 11L
  code[c3]           <- 12L
  code[f & p]        <- 13L
  code[f & c3]       <- 14L
  code[p & c3]       <- 15L
  code[f & p & c3]   <- 16L

  r <- terra::rast(depth)
  terra::values(r) <- code
  r[is.na(depth)] <- NA
  names(r) <- "gfisher_regions"

  if (verbose) print(table(terra::values(r, mat = FALSE), useNA = "ifany"))

  if (!is.null(dir.maps))
    fn.write_region_asc(r, file.path(dir.maps,
      sprintf("gfisher_survey_regions_%dmin.asc", res.min)))
  r
}


#' Combine the age-0 and GFISHER region grids.
#'
#' Age-0 regions (1-9) override GFISHER survey codes (10-16) where they overlap.
#'
#' @param age0 SpatRaster of age-0 region IDs, or a path to one. Pass the
#'   hand-edited `age0_survey_regions_5min_mod.asc` here if you have one --
#'   fn.make_age0_region_grids() regenerates the un-edited version.
#' @param gfisher SpatRaster of GFISHER survey codes, or a path to one.
#' @param dir.maps Output directory. NULL skips writing.
#' @param verbose Report the override count and the cell tally.
#' @return SpatRaster of combined region codes.
fn.combine_regions <- function(age0, gfisher, dir.maps = NULL, verbose = TRUE) {

  if (is.character(age0))    age0    <- terra::rast(age0)
  if (is.character(gfisher)) gfisher <- terra::rast(gfisher)
  if (!terra::compareGeom(age0, gfisher, stopOnError = FALSE))
    stop("age0 and gfisher grids do not share a geometry.")

  av <- terra::values(age0,    mat = FALSE)
  gv <- terra::values(gfisher, mat = FALSE)

  comb <- gv
  a0 <- !is.na(av) & av >= 1 & av <= 9
  comb[a0] <- av[a0]
  comb[is.na(av) | is.na(gv)] <- NA          # land stays NoData

  r <- terra::rast(gfisher)
  terra::values(r) <- comb
  names(r) <- "regions"

  if (verbose) {
    message("  age-0 cells overriding a survey code: ",
            sum(a0 & !is.na(gv) & gv >= 10))
    print(table(terra::values(r, mat = FALSE), useNA = "ifany"))
  }

  if (!is.null(dir.maps)) {
    res.min <- round(terra::res(r)[1] * 60, 0)
    fn.write_region_asc(r, file.path(dir.maps,
      sprintf("combined_regions_%dmin.asc", res.min)))
  }
  r
}


#' Write a region grid as an Ecospace ASCII, integer-coded.
#'
#' @param r SpatRaster of integer codes.
#' @param file Output path.
#' @return The path written, invisibly.
fn.write_region_asc <- function(r, file) {
  if (!dir.exists(dirname(file))) dir.create(dirname(file), recursive = TRUE)
  terra::writeRaster(r, file, overwrite = TRUE, NAflag = -9999,
                     datatype = "INT2S", filetype = "AAIGrid")
  unlink(paste0(file, ".aux.xml"))
  invisible(file)
}


#' Region code labels used by the attribute table and the plot legend.
REGION_LABELS <- c(
  "0" = "Water (unsampled)",
  "1" = "Saint Andrew Bay", "2" = "Saint Joe Bay",    "3" = "Turkey Point",
  "4" = "Mid Big Bend",     "5" = "Cedar Key",        "6" = "Tampa Bay",
  "7" = "Sarasota Bay",     "8" = "Charlotte Harbor", "9" = "Marco Island",
  "10" = "FWRI",            "11" = "PASC",            "12" = "PC",
  "13" = "FWRI+PASC",       "14" = "FWRI+PC",         "15" = "PASC+PC",
  "16" = "FWRI+PASC+PC")


#' Attribute table for a region grid.
#'
#' Areas are geodesic per-cell areas, so they account for a 5-arcmin cell
#' shrinking with latitude.
#'
#' @param r SpatRaster of region codes, or a path to one.
#' @param dir.maps Output directory. NULL skips writing.
#' @param include.land Add a -9999 "Land" row, as the legacy table had, so the
#'   cell counts sum to the full grid rather than to the water cells only.
#' @param verbose Print the table.
#' @return data.frame with ID, NAME, N_CELLS, AREA_KM2, AREA_ACRES.
fn.region_attributes <- function(r, dir.maps = NULL, include.land = TRUE,
                                 verbose = TRUE) {

  if (is.character(r)) r <- terra::rast(r)
  res.min <- round(terra::res(r)[1] * 60, 0)

  v   <- terra::values(r, mat = FALSE)
  akm <- terra::values(terra::cellSize(r, unit = "km"), mat = FALSE)

  vals <- sort(unique(v[!is.na(v)]))
  df <- data.frame(
    ID         = vals,
    NAME       = unname(REGION_LABELS[as.character(vals)]),
    N_CELLS    = vapply(vals, function(x) sum(v == x, na.rm = TRUE), integer(1)),
    AREA_KM2   = vapply(vals, function(x) round(sum(akm[which(v == x)]), 2), numeric(1)),
    stringsAsFactors = FALSE)

  # Land is NA in the raster but was a -9999 row in the legacy table.
  if (isTRUE(include.land) && any(is.na(v))) {
    land <- data.frame(ID = -9999L, NAME = "Land", N_CELLS = sum(is.na(v)),
                       AREA_KM2 = round(sum(akm[which(is.na(v))]), 2),
                       stringsAsFactors = FALSE)
    df <- rbind(land, df)
  }
  df$AREA_ACRES <- round(df$AREA_KM2 * 247.105381, 1)

  if (verbose) {
    print(df, row.names = FALSE)
    message(sprintf("  total cells: %d   total area: %.1f km2",
                    sum(df$N_CELLS), sum(df$AREA_KM2)))
  }

  if (!is.null(dir.maps)) {
    if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)
    utils::write.csv(df, file.path(dir.maps,
      sprintf("combined_regions_%dmin_attributes.csv", res.min)), row.names = FALSE)
  }
  df
}


#' Plot a combined region grid with a categorical legend.
#'
#' @param r SpatRaster of region codes, or a path to one.
#' @param dir.maps Output directory for the PNG.
#' @param tag Optional filename suffix.
#' @return The path written, invisibly.
fn.plot_regions <- function(r, dir.maps, tag = "") {

  if (is.character(r)) r <- terra::rast(r)
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)
  res.min <- round(terra::res(r)[1] * 60, 0)

  cols <- c("#e8f2fb",                                          # 0  water
            "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",  # 1-5  age0
            "#8c564b","#e377c2","#7f7f7f","#bcbd22",            # 6-9  age0
            "#17becf","#aec7e8","#ffbb78","#98df8a","#ff9896",  # 10-14 surveys
            "#c5b0d5","#c49c94")                                # 15-16 surveys
  names(cols) <- as.character(0:16)

  present <- sort(unique(terra::values(r, mat = FALSE)))
  present <- present[!is.na(present)]
  labs <- paste(present, unname(REGION_LABELS[as.character(present)]))

  fout <- file.path(dir.maps, sprintf("combined_regions_%dmin%s.png", res.min,
                                      if (nzchar(tag)) paste0(" ", tag) else ""))
  png(fout, width = 10, height = 9, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)

  terra::plot(r, type = "classes", levels = present,
              col = unname(cols[as.character(present)]),
              colNA = "#d9cdb8", mar = c(4, 4, 3, 12),
              plg = list(legend = labs, cex = 0.8),
              main = sprintf("Combined survey regions (%d min)", res.min))
  maps::map("state", add = TRUE, col = "gray30", lwd = 0.5)
  box()

  message("Figure written to\n  ", fout)
  invisible(fout)
}
