#' Fisheries management area basemaps for WFS Ecospace.
#'
#' Rewrite of "maps/management areas/make management area maps.R" as callable
#' functions. Each zipped shapefile becomes one Ecospace ASCII grid aligned to
#' the depth template:
#'
#'   1       cells the management area covers
#'   0       other cells inside the model domain
#'   -9999   land and excluded deep water (NODATA in the depth grid)
#'
#' One zip can yield more than one area: `madswan_steamboat_edges` holds three
#' features in a LABEL field, and Madison/Swanson and Steamboat Lumps share
#' seasonal management while the Edges differ, so it is split into two grids.
#' See `fn.default_ma_splits()`.
#'
#' terra only. Do not add library(raster) or library(sp) here -- see the note at
#' the top of R/GFISHER functions.R.


#' Default split rules for zips holding more than one management area.
#'
#' Each entry names a zip (without extension) and gives the attribute field to
#' group by, plus the output name for each group of field values.
#'
#' @return Named list of split specifications.
fn.default_ma_splits <- function() {
  list(
    madswan_steamboat_edges = list(
      field  = "LABEL",
      groups = list(
        "madison_swanson_steamboat" = c("Madison and Swanson sites",
                                        "Steamboat Lumps"),
        "edges"                     = c("Edges")
      )
    )
  )
}


#' Pick the most useful shapefile from an unzipped folder.
#'
#' Prefers polygon (_po), then line (_ln), then point (_pt); otherwise the first
#' shapefile found. `.shp.xml` sidecars are ignored.
#'
#' @param dir Directory to search, recursively.
#' @return Path to one .shp, or NA_character_ if none.
fn.pick_shapefile <- function(dir) {
  shps <- list.files(dir, pattern = "\\.shp$", full.names = TRUE,
                     recursive = TRUE, ignore.case = TRUE)
  shps <- shps[!grepl("\\.shp\\.xml$", shps, ignore.case = TRUE)]
  if (length(shps) == 0) return(NA_character_)
  for (tag in c("_po", "_ln", "_pt")) {
    hit <- shps[grepl(tag, basename(shps), ignore.case = TRUE)]
    if (length(hit)) return(hit[1])
  }
  shps[1]
}


#' Rasterize one management area onto the depth grid.
#'
#' @param v SpatVector of the area, already in the depth CRS.
#' @param depth SpatRaster template.
#' @param inside.val,outside.val Values for covered and uncovered model cells.
#' @param touches Include any cell the area touches, not only those whose centre
#'   it contains. TRUE is the legacy behaviour and is the safer choice for thin
#'   or small areas that would otherwise vanish at coarse resolution.
#' @return Single-layer SpatRaster.
fn.rasterize_management_area <- function(v, depth, inside.val = 1,
                                         outside.val = 0, touches = TRUE) {

  r <- terra::rasterize(v, depth, field = 1, background = NA, touches = touches)

  out <- terra::rast(depth)
  terra::values(out) <- outside.val
  out[!is.na(r)]     <- inside.val
  out[is.na(depth)]  <- NA
  out
}


#' Build management area basemaps from a folder of zipped shapefiles.
#'
#' Entry point called from make_WFS_basemaps.R.
#'
#' @param depth SpatRaster of positive depths (the template grid).
#' @param dir.zips Directory holding the zipped shapefiles.
#' @param dir.maps Output directory for ASCII and plots.
#' @param splits Split rules; see fn.default_ma_splits(). NULL for none.
#' @param inside.val,outside.val Values for covered and uncovered model cells.
#' @param nodata.val NODATA flag written to the ASCII header.
#' @param touches Passed to terra::rasterize().
#' @param write.ascii Write Ecospace ASCII grids.
#' @param verbose Report each area and its cell count.
#' @return SpatRaster, one layer per management area.
fn.make_management_area_maps <- function(depth, dir.zips, dir.maps,
                                         splits      = fn.default_ma_splits(),
                                         inside.val  = 1,
                                         outside.val = 0,
                                         nodata.val  = -9999,
                                         touches     = TRUE,
                                         write.ascii = TRUE,
                                         verbose     = TRUE) {

  if (inherits(depth, "Raster")) depth <- terra::rast(depth)
  stopifnot(inherits(depth, "SpatRaster"))
  # The depth ascii is lon/lat WGS84 but carries no .prj, so stamp it before
  # anything gets projected onto it.
  if (is.na(terra::crs(depth)) || !nzchar(terra::crs(depth)))
    terra::crs(depth) <- "EPSG:4326"

  res.min <- round(terra::res(depth)[1] * 60, 0)
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)

  zips <- list.files(dir.zips, pattern = "\\.zip$", full.names = TRUE)
  if (length(zips) == 0) stop("No .zip files found in ", dir.zips)
  if (verbose) message("Found ", length(zips), " zipped shapefile(s).")

  lyrs <- list()
  for (z in zips) {

    area <- tools::file_path_sans_ext(basename(z))
    if (verbose) message("---- ", area, " ----")

    ex.dir <- file.path(tempdir(), paste0("ma_", area))
    dir.create(ex.dir, showWarnings = FALSE, recursive = TRUE)
    utils::unzip(z, exdir = ex.dir, overwrite = TRUE)

    shp <- fn.pick_shapefile(ex.dir)
    if (is.na(shp)) {
      warning("No shapefile inside ", basename(z), "; skipping.")
      next
    }
    if (verbose) message("  using: ", basename(shp))

    v <- terra::vect(shp)
    if (is.na(terra::crs(v)) || !nzchar(terra::crs(v)))
      terra::crs(v) <- "EPSG:4326"
    v <- terra::project(v, depth)

    # One zip, one area -- unless a split rule says otherwise.
    todo <- list()
    if (!is.null(splits) && area %in% names(splits)) {
      spec <- splits[[area]]
      if (!spec$field %in% names(v))
        stop("Split field '", spec$field, "' not found in ", basename(shp),
             ". Fields: ", paste(names(v), collapse = ", "))
      for (g in names(spec$groups)) {
        vg <- v[v[[spec$field]][, 1] %in% spec$groups[[g]], ]
        if (verbose) message(sprintf("  group '%s' (%d features)", g, nrow(vg)))
        if (nrow(vg) == 0) {
          warning("No features for group '", g, "' in ", area, "; skipping.")
          next
        }
        todo[[g]] <- vg
      }
    } else {
      todo[[area]] <- v
    }

    for (nm in names(todo)) {
      r <- fn.rasterize_management_area(todo[[nm]], depth, inside.val,
                                        outside.val, touches)
      names(r) <- nm
      n.in <- sum(terra::values(r, mat = FALSE) == inside.val, na.rm = TRUE)
      if (verbose) message("    cells inside area: ", n.in)
      if (n.in == 0)
        warning("Area '", nm, "' covers no model cells at ", res.min,
                " min; it may fall outside the grid or be too small to touch a cell.")
      lyrs[[nm]] <- r
    }
  }

  if (length(lyrs) == 0) stop("No management areas were rasterized.")
  ma <- terra::rast(lyrs)
  names(ma) <- names(lyrs)

  if (isTRUE(write.ascii)) {
    for (nm in names(ma)) {
      f <- file.path(dir.maps, sprintf("%s_%dmin.asc", nm, res.min))
      terra::writeRaster(ma[[nm]], f, overwrite = TRUE, NAflag = nodata.val,
                         datatype = "INT2S", filetype = "AAIGrid")
      unlink(paste0(f, ".aux.xml"))
    }
    if (verbose) message("Ecospace ascii files written to\n  ", dir.maps)
  }

  invisible(ma)
}


#' Plot the management area grids as a panel.
#'
#' Binary layers, so a two-colour scale rather than the continuous ramp used for
#' the habitat maps.
#'
#' @param x SpatRaster from fn.make_management_area_maps(), or a directory of
#'   written ASCII grids.
#' @param dir.maps Output directory for the PNG.
#' @param tag Optional filename suffix.
#' @param col Length-2 colour vector: outside, inside.
#' @return The path written, invisibly.
fn.plot_management_areas <- function(x, dir.maps, tag = "",
                                     col = c("gray95", "firebrick")) {

  if (is.character(x)) {
    f <- list.files(x, pattern = "\\.asc$", full.names = TRUE)
    x <- terra::rast(f)
    names(x) <- sub("_[0-9]+min$", "", tools::file_path_sans_ext(basename(f)))
  }
  if (missing(dir.maps)) stop("dir.maps is required.")
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)

  res.min <- round(terra::res(x)[1] * 60, 0)
  n  <- terra::nlyr(x)
  nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)

  fout <- file.path(dir.maps, sprintf("management areas %dmin%s.png", res.min,
                                      if (nzchar(tag)) paste0(" ", tag) else ""))
  png(fout, height = 3.2 * nr, width = 3.6 * nc, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(nr, nc))

  for (i in seq_len(n)) {
    v    <- terra::values(x[[i]], mat = FALSE)
    n.in <- sum(v == 1, na.rm = TRUE)

    # A layer with no cells inside is constant, and terra drops `breaks` for a
    # constant raster and uses the last colour -- which painted empty areas
    # solid red. Give those a single explicit colour instead.
    if (n.in == 0) {
      terra::plot(x[[i]], col = col[1], colNA = "lightgray", legend = FALSE,
                  mar = c(2, 2, 3.5, 1),
                  main = sprintf("%s\n0 cells (outside the grid)", names(x)[i]))
    } else {
      terra::plot(x[[i]], col = col, breaks = c(-0.5, 0.5, 1.5),
                  colNA = "lightgray", legend = FALSE,
                  mar = c(2, 2, 3.5, 1),
                  main = sprintf("%s\n%d cells", names(x)[i], n.in))
    }
    maps::map("state", add = TRUE, fill = TRUE, col = "lightgray")
  }

  message("Figure written to\n  ", fout)
  invisible(fout)
}
