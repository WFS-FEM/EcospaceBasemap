

#' Download the Gulf Data Atlas seagrass layer.
#'
#' Fetches GulfwideSAV.zip from NOAA NCEI and unpacks it into a subdirectory
#' named after the archive. The FWC Seagrass_Statewide layer is not downloadable
#' this way and has to be supplied by hand.
#'
#' @param dir.out Directory to download into; a GulfwideSAV/ subdirectory is
#'   created inside it.
#' @return NULL, invisibly. Called for the download.
fn.pull_seagrass <- function(dir.out){
  
  #dir.out = "C:\\Users\\dchagaris\\Github\\WFS-FEM\\EnvironmentalDrivers2EwE\\data\\seagrass"
  
  message('Downloading seagrass layer from Gulf Data Atlas\nhttps://www.ncei.noaa.gov/maps/gulf-data-atlas/atlas.htm') 
  
  url = "https://www.ncei.noaa.gov/waf/data-atlas-waf/biotic/documents/GulfwideSAV.zip"
  file.out = file.path(dir.out,basename(url))
  download.file(url, destfile = file.out, mode='wb')
  unzip(zipfile=file.out, exdir=file.path(dir.out, gsub(".zip","",basename(url))))
  unlink(file.out)

  message('Seagrass data downloaded and extracted to \n',dir.out)
} #eof

#' Rasterize the union of two or more seagrass polygon sets.
#'
#' The Gulf Data Atlas (GulfwideSAV) and FWC (Seagrass_Statewide) layers map the
#' same beds, not different ones -- at 5 min, 186 of 247 Gulfwide cells and 186
#' of 212 FWC cells overlap, and summing the finished rasters would push 139
#' cells past 1. Combining has to start from the polygons, where an overlap can
#' be resolved, rather than from the grids, where it cannot be distinguished
#' from two adjacent beds in one cell.
#'
#' Note the two sources differ in vintage: GulfwideSAV carries `hab_88`/`hab_92`
#' fields and reports roughly twice the coverage of the newer FWC layer where
#' both are present. The union is therefore closer to historical maximum extent
#' than to current extent. Pass a single directory to get one source alone.
#'
#' @param dirs.seagrass Character vector of directories, each holding one .shp.
#' @param dir.ascii Output directory. NULL skips writing.
#' @param depth SpatRaster template.
#' @param label Name used in the output filenames.
#' @param make.valid Repair invalid geometries first. GulfwideSAV reports a
#'   whole-globe bounding box, so at least one feature is corrupt.
#' @param fact Sub-grid factor used to estimate coverage. Each model cell is
#'   built from fact^2 sub-cells; 20 gives roughly 460 m sub-cells at 5 min.
#'   Larger is more accurate and slower (cost is O(fact^2)).
#' @param verbose Report feature counts through each step.
#' @return SpatRaster of the fraction of each cell covered by any source.
fn.combine_seagrass <- function(dirs.seagrass, dir.ascii, depth,
                                label = "union", make.valid = TRUE,
                                fact = 20, verbose = TRUE) {

  if (inherits(depth, "Raster")) depth <- terra::rast(depth)
  if (is.na(terra::crs(depth)) || !nzchar(terra::crs(depth)))
    terra::crs(depth) <- "EPSG:4326"
  res.min <- round(terra::res(depth)[1] * 60, 0)

  parts <- list()
  for (d in dirs.seagrass) {
    f <- list.files(d, pattern = "\\.shp$", full.names = TRUE)
    if (length(f) == 0) stop("No .shp found in ", d)
    v <- terra::vect(f[1])
    if (verbose) message(sprintf("  %-20s %6d features", basename(d), nrow(v)))

    if (isTRUE(make.valid)) {
      bad <- !terra::is.valid(v)
      if (any(bad)) {
        if (verbose) message(sprintf("    repairing %d invalid geometr%s",
                                     sum(bad), if (sum(bad) == 1) "y" else "ies"))
        v <- terra::makeValid(v)
      }
    }

    v <- terra::project(v, terra::crs(depth))
    # Crop early: it drops everything outside the model domain and takes the
    # corrupt whole-globe features with it.
    v <- terra::crop(v, terra::ext(depth))
    if (verbose) message(sprintf("    %6d features inside the grid", nrow(v)))

    if (nrow(v) > 0) parts[[length(parts) + 1]] <- v[, 0]   # geometry only
  }
  if (length(parts) == 0) stop("No seagrass features inside the model domain.")

  # Union by presence on a sub-grid, rather than by dissolving polygons.
  #
  # A fine cell is either inside some polygon or it is not, so presence is
  # inherently a union -- overlapping beds cannot be counted twice, and the
  # result cannot exceed 1 by construction. Averaging presence back up to the
  # model grid gives the covered fraction.
  #
  # The alternatives both failed on this data. terra::aggregate() groups
  # geometries without merging overlapping rings, so rasterize(cover = TRUE)
  # counted shared area once per ring and pushed 251 of 273 cells above 1.
  # sf::st_union() on the ~74k repaired polygons took 15 minutes and still
  # produced overlaps. This runs in seconds.
  #
  # Accuracy is set by `fact`: each model cell is estimated from fact^2 sub-cells
  # classified by whether their centre falls inside a bed.
  v <- if (length(parts) == 1) parts[[1]] else do.call(rbind, parts)

  if (verbose) message(sprintf("  rasterizing presence on a %dx sub-grid...", fact))
  fine <- terra::disagg(terra::rast(depth), fact = fact)
  pres <- terra::rasterize(v, fine, field = 1, background = 0)
  sav  <- terra::aggregate(pres, fact = fact, fun = "mean", na.rm = TRUE)
  sav  <- terra::resample(sav, depth, method = "near")   # guard exact alignment
  sav[is.na(sav)] <- 0
  sav[is.na(depth)] <- NA
  names(sav) <- "seagrass"

  if (verbose) {
    v <- terra::values(sav, mat = FALSE)
    message(sprintf("  cells > 0: %d   mean %.5f   max %.4f",
                    sum(v > 0, na.rm = TRUE), mean(v, na.rm = TRUE),
                    max(v, na.rm = TRUE)))
  }

  if (!is.null(dir.ascii)) {
    if (!dir.exists(dir.ascii)) dir.create(dir.ascii, recursive = TRUE)
    f <- file.path(dir.ascii, sprintf("seagrass_coverage_%s_%dmin.asc", label, res.min))
    terra::writeRaster(sav, f, overwrite = TRUE, NAflag = -9999,
                       gdal = c("DECIMAL_PRECISION=4"))
    unlink(paste0(f, ".aux.xml"))

    png(gsub("\\.asc$", ".png", f), height = 7, width = 7, units = "in", res = 300)
    terra::plot(sav, colNA = "lightgray", mar = c(3, 3, 3, 6),
                col = colorRampPalette(colorRamps::matlab.like2(100),
                                       bias = 2, interpolate = "spline")(100),
                main = paste0("Seagrass proportion coverage\nunion of: ",
                              paste(basename(dirs.seagrass), collapse = " + ")))
    maps::map("state", add = TRUE, fill = TRUE, col = "lightgray")
    dev.off()
    message("Seagrass union written to\n  ", dir.ascii)
  }

  sav
}


#' Combine seagrass coverage rasters without going back to the polygons.
#'
#' A shortcut for when the per-source grids are already in the workspace.
#'
#' The union cannot be recovered exactly from coverage rasters. For a cell
#' holding `a` and `b`, the truth is `a + b - overlap`, and the overlap is not
#' recorded, so all that is known is:
#'
#'   max(a, b)  <=  union  <=  min(1, a + b)
#'
#' `"max"` returns the lower bound, `"sum_capped"` the upper. Running both
#' brackets the answer and shows how much the choice matters.
#'
#' For the WFS pair, `"max"` should be close to exact: the two sources map the
#' same beds, with FWC nested inside Gulfwide in 168 of the 186 cells where both
#' are non-zero (means 0.824 vs 0.386). When one source is a near-subset of the
#' other the overlap approaches `min(a, b)`, and `a + b - min(a, b)` is exactly
#' `max(a, b)`. Worth confirming against fn.combine_seagrass() once; if the two
#' agree, this is the cheaper call to keep using, since the cost of the polygon
#' route is reading and repairing ~90k features, not the rasterizing.
#'
#' @param x Multi-layer SpatRaster, or a list of single-layer SpatRasters, of
#'   coverage proportions on a common grid.
#' @param method "max" for the lower bound, "sum_capped" for the upper.
#' @param dir.ascii Output directory. NULL skips writing.
#' @param depth Optional template; used to enforce the land mask.
#' @param label Name used in the output filenames. Defaults to `method`.
#' @param verbose Report the per-source and combined totals.
#' @return Single-layer SpatRaster named "seagrass".
fn.combine_seagrass_rasters <- function(x, method = c("max", "sum_capped"),
                                        dir.ascii = NULL, depth = NULL,
                                        label = NULL, verbose = TRUE) {

  method <- match.arg(method)
  if (is.list(x)) x <- terra::rast(x)
  stopifnot(inherits(x, "SpatRaster"))
  if (terra::nlyr(x) < 2)
    stop("Need at least two layers to combine; got ", terra::nlyr(x), ".")
  if (is.null(label)) label <- method

  sav <- switch(method,
                max        = max(x, na.rm = TRUE),
                sum_capped = min(sum(x, na.rm = TRUE), 1))

  if (!is.null(depth)) {
    if (inherits(depth, "Raster")) depth <- terra::rast(depth)
    sav[is.na(sav) & !is.na(depth)] <- 0
    sav[is.na(depth)] <- NA
  }
  names(sav) <- "seagrass"

  if (verbose) {
    for (i in seq_len(terra::nlyr(x))) {
      v <- terra::values(x[[i]], mat = FALSE)
      message(sprintf("  %-24s n>0 %4d  mean %.5f  total %.2f",
                      names(x)[i], sum(v > 0, na.rm = TRUE),
                      mean(v, na.rm = TRUE), sum(v, na.rm = TRUE)))
    }
    v <- terra::values(sav, mat = FALSE)
    message(sprintf("  %-24s n>0 %4d  mean %.5f  total %.2f",
                    paste0("combined (", method, ")"), sum(v > 0, na.rm = TRUE),
                    mean(v, na.rm = TRUE), sum(v, na.rm = TRUE)))
    n.over <- sum(v > 1 + 1e-9, na.rm = TRUE)
    if (n.over > 0) warning(n.over, " cell(s) exceed 1.")
  }

  if (!is.null(dir.ascii)) {
    if (!dir.exists(dir.ascii)) dir.create(dir.ascii, recursive = TRUE)
    res.min <- round(terra::res(sav)[1] * 60, 0)
    f <- file.path(dir.ascii, sprintf("seagrass_coverage_%s_%dmin.asc", label, res.min))
    terra::writeRaster(sav, f, overwrite = TRUE, NAflag = -9999,
                       gdal = c("DECIMAL_PRECISION=4"))
    unlink(paste0(f, ".aux.xml"))

    png(gsub("\\.asc$", ".png", f), height = 7, width = 7, units = "in", res = 300)
    terra::plot(sav, colNA = "lightgray", mar = c(3, 3, 3, 6),
                col = colorRampPalette(colorRamps::matlab.like2(100),
                                       bias = 2, interpolate = "spline")(100),
                main = sprintf("Seagrass proportion coverage\ncombined (%s) of %d sources",
                               method, terra::nlyr(x)))
    maps::map("state", add = TRUE, fill = TRUE, col = "lightgray")
    dev.off()

    message("Combined seagrass written to\n  ", f)
  }

  sav
}


#' Rasterize ONE seagrass polygon source onto the Ecospace grid.
#'
#' `cover = TRUE` returns the fraction of each cell covered by polygons, so the
#' output is a proportion in [0, 1]. Cells with no polygon become 0 and land
#' becomes NA.
#'
#' Use this per source to inspect the layers individually. To get a single map
#' for the sum-to-1 step, combine the sources with fn.combine_seagrass() or
#' fn.combine_seagrass_rasters() -- the two available sources map the same beds,
#' so adding their rasters double-counts.
#'
#' @param dir.seagrass Directory holding one .shp.
#' @param dir.ascii Output directory for the ASCII grid and PNG.
#' @param depth SpatRaster template (the model grid).
#' @return Single-layer SpatRaster of proportional cover.
fn.make_seagrass_ascii <- function(dir.seagrass, dir.ascii, depth){
  
  #dir.seagrass= "C:\\Users\\dchagaris\\Github\\WFS-FEM\\EnvironmentalDrivers2EwE\\data\\seagrass\\Seagrass_Statewide" 
  #dir.ascii = file.path(dir.maps,"seagrass")
  #depth = depth.15min
  
  res = res(depth)[1]*60

  # Load required packages
  require(terra)
  
  # Read the shapefile
  file.shp <- list.files(dir.seagrass, pattern=".shp$", full.names=T)
  seagrass.shp <- vect(file.shp)
  #seagrass.shp <- vect(file.path(dir.seagrass, "Seagrass_ALFLMSTX.shp"))

  # Ensure same CRS as depth raster
  seagrass.shp <- project(seagrass.shp, crs(depth))
  

  #plot(seagrass.shp)
  # Rasterize with weights=TRUE calculates the fraction of each cell covered by polygons
  seagrass.prop <- rasterize(seagrass.shp, depth, cover=TRUE, touches=T)
  seagrass.prop[is.na(seagrass.prop)] <- 0
  seagrass.prop[is.na(depth)] <- NA
  plot(seagrass.prop, col=c('white',matlab.like(100)), colNA='gray')
  # Create output directory if it doesn't exist
  if (!dir.exists(dir.ascii)) dir.create(dir.ascii, recursive=TRUE)
  
  # Write output to ASCII file
  terra::writeRaster(seagrass.prop, filename=file.path(dir.ascii, paste0("seagrass_coverage_",basename(dir.seagrass),'_',res,"min.asc")),
                     overwrite=T, gdal=c('DECIMAL_PRECISION=4'), NAflag=-9999)
  #plot
  datasource = ifelse(basename(dir.seagrass)=='Seagrass_Statewide','FWC','Gulf Data Atlas')
  png(file.path(dir.ascii, paste0("seagrass_coverage_",basename(dir.seagrass),'_',res,"min.png")),height = 7, width=7, units='in', res=300)
  plot(seagrass.prop, colNA='gray', main=paste0('Seagrass Proportion Coverage\nsource: ',datasource))
  dev.off()
  
  return(seagrass.prop)
}#eof
