#' Fleet port basemaps for WFS Ecospace.
#'
#' Consolidates the eight scripts under "maps/ports" into one parameterised
#' pipeline. Every fleet follows the same skeleton, differing only in its data
#' source, species filter, majority threshold and value column:
#'
#'   1. rank counties by landings / trips / vessels
#'   2. take the smallest set of counties reaching `cum.thresh` percent
#'   3. give each of those counties one port cell
#'   4. write a binary grid: 1 = port, 0 = water, -9999 = all other land
#'
#' A port is a COASTAL LAND cell of the depth grid:
#'   land     - NoData in the depth grid, which is what defines land here
#'   coastal  - touches at least one water cell (queen adjacency)
#'   in-county- cell centre inside the county polygon
#'   nearest  - of those, closest to the county's anchor point (see POP_CENTER)
#'
#' Gulf vs Atlantic is decided by an explicit county list, not by distance,
#' because the grid's eastern edge clips Atlantic water near Jacksonville and
#' Cape Canaveral.
#'
#' terra + sf only. The legacy coastal_ports.R was raster/sp-based; see the note
#' at the top of R/GFISHER functions.R for why that matters here.


# ---- reference tables -------------------------------------------------------

#' Gulf of Mexico / West Florida coastal counties covered by the grid.
GULF_COUNTIES <- c("ESCAMBIA","SANTA ROSA","OKALOOSA","WALTON","BAY","GULF",
  "FRANKLIN","WAKULLA","JEFFERSON","TAYLOR","DIXIE","LEVY","CITRUS","HERNANDO",
  "PASCO","PINELLAS","HILLSBOROUGH","MANATEE","SARASOTA","CHARLOTTE","LEE",
  "COLLIER","MONROE")

#' Largest population centre of each Gulf county (WGS84).
#'
#' The port cell is the coastal land cell in that county nearest this point, so
#' an inland county seat still resolves to the right stretch of coast --
#' Crestview to Fort Walton Beach, Chiefland to Cedar Key, Perry to the Big Bend.
#' Edit this table to move a port.
POP_CENTER <- utils::read.csv(text =
"County,Town,lon,lat
ESCAMBIA,Pensacola,-87.2169,30.4213
SANTA ROSA,Navarre,-86.8635,30.4013
OKALOOSA,Crestview,-86.5705,30.7621
WALTON,DeFuniak Springs,-86.1150,30.7213
BAY,Panama City,-85.6602,30.1588
GULF,Port St. Joe,-85.3018,29.8135
FRANKLIN,Apalachicola,-84.9856,29.7258
WAKULLA,Crawfordville,-84.3752,30.1758
JEFFERSON,Monticello,-83.8710,30.5455
TAYLOR,Perry,-83.5827,30.1174
DIXIE,Cross City,-83.1257,29.6352
LEVY,Chiefland,-82.8598,29.4747
CITRUS,Homosassa Springs,-82.5757,28.8003
HERNANDO,Spring Hill,-82.5254,28.4769
PASCO,New Port Richey,-82.7193,28.2442
PINELLAS,St. Petersburg,-82.6403,27.7676
HILLSBOROUGH,Tampa,-82.4572,27.9506
MANATEE,Bradenton,-82.5748,27.4989
SARASOTA,Sarasota,-82.5307,27.3364
CHARLOTTE,Port Charlotte,-82.0906,26.9762
LEE,Cape Coral,-81.9495,26.5629
COLLIER,Naples,-81.7948,26.1420
MONROE,Key West,-81.7800,24.5551
", stringsAsFactors = FALSE)

#' Headboat vessels by county (scoped 2026-08-19).
#'
#' Monroe is excluded: the Keys headboats fish the Atlantic reef, outside this
#' shelf grid. Lee is counted at its pre-Hurricane-Ian level; set to 0 for a
#' present-day fleet.
HEADBOATS <- utils::read.csv(text =
"County,Boats
OKALOOSA,10
PINELLAS,6
BAY,3
MANATEE,2
SARASOTA,2
COLLIER,2
LEE,2
ESCAMBIA,1
FRANKLIN,1
GULF,0
WALTON,0
WAKULLA,0
TAYLOR,0
DIXIE,0
LEVY,0
CITRUS,0
HERNANDO,0
PASCO,0
HILLSBOROUGH,0
CHARLOTTE,0
", stringsAsFactors = FALSE)

#' Headboat home marinas, which override POP_CENTER for that fleet.
#'
#' Headboats are tied to a known dock rather than to a county population centre.
#' Counties absent here fall back to POP_CENTER.
HEADBOAT_HUB <- utils::read.csv(text =
"County,Town,lon,lat
OKALOOSA,Destin Harbor,-86.5133,30.3935
PINELLAS,Clearwater Beach Marina,-82.8270,27.9764
BAY,Capt Anderson's Marina PCB,-85.7530,30.1466
MANATEE,Cortez fishing village,-82.6890,27.4677
SARASOTA,Marina Jack Sarasota,-82.5477,27.3345
COLLIER,Naples City Dock,-81.7920,26.1387
LEE,Getaway Marina Ft Myers Bch,-81.9502,26.4548
ESCAMBIA,Perdido Key,-87.4500,30.3080
FRANKLIN,Apalachicola waterfront,-84.9820,29.7250
", stringsAsFactors = FALSE)


# ---- data readers -----------------------------------------------------------

#' Read and clean the FWC commercial landings summary.
#'
#' @param file Path to ReportCreatorResults-County.csv (9 header lines).
#' @return data.frame with Year, County, Species, Pounds, Trips, AvgPrice, Value.
fn.read_fwc_landings <- function(file) {
  d <- utils::read.csv(file, skip = 9, stringsAsFactors = FALSE,
                       check.names = FALSE)
  names(d) <- c("Year", "County", "Species", "Pounds", "Trips", "AvgPrice", "Value")
  d$Pounds <- suppressWarnings(as.numeric(gsub(",", "", d$Pounds)))
  d <- d[!is.na(d$Year) & d$County != "" & !is.na(d$Pounds), ]
  d$County <- toupper(trimws(d$County))
  # Not a mappable county.
  d[d$County != "INLAND/OUT-OF-STATE", ]
}


#' Read MRIP directed trips and resolve FIPS codes to county names.
#'
#' @param file Path to the MRIP dtrips CSV.
#' @param fl.fips Florida state FIPS; county codes are within-state.
#' @return data.frame with County, common, mode_fx, dtrip.
fn.read_mrip_dtrips <- function(file, fl.fips = 12) {
  m <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  m$dtrip  <- suppressWarnings(as.numeric(m$dtrip))
  m <- m[!is.na(m$dtrip) & !is.na(m$cnty), ]
  m$common <- toupper(trimws(m$common))

  utils::data("county.fips", package = "maps", envir = environment())
  cf <- get("county.fips", envir = environment())
  fips <- cf[grepl("^florida,", cf$polyname), ]
  fips$County <- toupper(sub(":.*", "", sub("florida,", "", fips$polyname)))
  fips$cnty   <- fips$fips - fl.fips * 1000
  fips <- unique(fips[, c("cnty", "County")])

  m$County <- fips$County[match(m$cnty, fips$cnty)]
  if (any(is.na(m$County)))
    warning("Unmatched MRIP FIPS codes: ",
            paste(unique(m$cnty[is.na(m$County)]), collapse = ", "))
  m[!is.na(m$County), ]
}


#' Florida county polygons as sf, with an uppercase County column.
fn.florida_counties <- function() {
  flc <- sf::st_as_sf(maps::map("county", "florida", plot = FALSE, fill = TRUE))
  # maps returns unprojected lon/lat but st_as_sf stamps its own CRS, which does
  # not compare equal to EPSG:4326 and makes st_intersects() fail later. Force
  # it. sf warns that this does not reproject; that is exactly the intent, so
  # the warning is suppressed rather than worked around.
  suppressWarnings(sf::st_crs(flc) <- 4326)
  flc$County <- toupper(sub("florida,", "", flc$ID))
  flc
}


# ---- port cell machinery ----------------------------------------------------

#' Find the coastal land cells of the depth grid, indexed by county.
#'
#' Land is NoData in the depth grid; a cell is coastal if it touches at least
#' one water cell under queen adjacency.
#'
#' @param depth SpatRaster template.
#' @param flc sf of Florida counties with a County column.
#' @param counties Counties to index.
#' @param verbose Report the cell tally.
#' @return list(depth, dv, cpts, cty_cells, coast_cell).
fn.build_coastal <- function(depth, flc, counties = GULF_COUNTIES,
                             verbose = TRUE) {

  dv <- terra::values(depth, mat = FALSE)
  land_cell <- which(is.na(dv))

  # terra::adjacent has no `target` argument, so pairs are filtered to those
  # whose neighbour is water.
  adj <- terra::adjacent(depth, cells = land_cell, directions = 8, pairs = TRUE)
  adj <- adj[!is.na(dv[adj[, 2]]), , drop = FALSE]

  coast_cell <- sort(unique(adj[, 1]))
  nwater <- table(factor(adj[, 1], levels = coast_cell))
  cxy <- terra::xyFromCell(depth, coast_cell)

  cpts <- sf::st_as_sf(data.frame(cell = coast_cell, nwater = as.integer(nwater),
                                  x = cxy[, 1], y = cxy[, 2]),
                       coords = c("x", "y"), crs = 4326)

  cty_cells <- stats::setNames(vector("list", length(counties)), counties)
  for (cty in counties) {
    poly <- flc[flc$County == cty, ]
    if (nrow(poly) == 0) next
    cty_cells[[cty]] <- which(lengths(sf::st_intersects(cpts, sf::st_union(poly))) > 0)
  }

  if (verbose)
    message(sprintf("  depth grid: %d water, %d land, %d coastal-land candidates",
                    sum(!is.na(dv)), length(land_cell), length(coast_cell)))

  list(depth = depth, dv = dv, cpts = cpts, cty_cells = cty_cells,
       coast_cell = coast_cell)
}


#' Blank port row, for counties with no in-grid coastal cell.
PORT_BLANK <- data.frame(Town = NA_character_, town_km = NA_real_,
                         cell = NA_integer_, nwater = NA_integer_,
                         row = NA_integer_, col = NA_integer_,
                         lon = NA_real_, lat = NA_real_)


#' Assign one county its port cell.
#'
#' @param cty County name.
#' @param CL Output of fn.build_coastal().
#' @param pop.center Anchor table; POP_CENTER, or a fleet-specific override.
#' @return One-row data.frame, or NULL if the county has no candidate cell.
fn.assign_port <- function(cty, CL, pop.center = POP_CENTER) {
  idx <- CL$cty_cells[[cty]]
  pc  <- pop.center[pop.center$County == cty, ]
  if (is.null(idx) || length(idx) == 0 || nrow(pc) == 0) return(NULL)

  town <- sf::st_sfc(sf::st_point(c(pc$lon, pc$lat)), crs = 4326)
  dd   <- as.numeric(sf::st_distance(CL$cpts[idx, ], town))
  j    <- idx[which.min(dd)]
  ll   <- sf::st_coordinates(CL$cpts[j, ])

  data.frame(Town = pc$Town, town_km = round(min(dd) / 1000, 2),
             cell = CL$cpts$cell[j], nwater = CL$cpts$nwater[j],
             row = terra::rowFromCell(CL$depth, CL$cpts$cell[j]),
             col = terra::colFromCell(CL$depth, CL$cpts$cell[j]),
             lon = round(ll[1], 4), lat = round(ll[2], 4),
             stringsAsFactors = FALSE)
}


#' Rank counties and flag the majority set.
#'
#' @param bc data.frame with County and the value column.
#' @param value.col Name of the value column.
#' @param cum.thresh Cumulative percent defining the majority set.
#' @return `bc` ordered descending, with pct, cum and majority added.
fn.county_majority <- function(bc, value.col = "Pounds", cum.thresh = 80) {
  bc <- bc[order(-bc[[value.col]]), ]
  bc$pct <- 100 * bc[[value.col]] / sum(bc[[value.col]])
  bc$cum <- cumsum(bc$pct)
  # Smallest set reaching the threshold: through the first county with cum >= it.
  n.major <- which(bc$cum >= cum.thresh)[1]
  if (is.na(n.major)) n.major <- nrow(bc)
  bc$majority <- seq_len(nrow(bc)) <= n.major
  attr(bc, "n_major") <- n.major
  bc
}


#' Build the port assignment table for the majority counties.
#'
#' @param maj Majority subset of `bc`.
#' @param bc Full ranked county table.
#' @param CL Output of fn.build_coastal().
#' @param value.col Name of the value column.
#' @param pop.center Anchor table.
#' @return One row per majority county, ordered by rank.
fn.port_table <- function(maj, bc, CL, value.col = "Pounds",
                          pop.center = POP_CENTER) {
  res <- data.frame()
  for (cty in maj$County) {
    in.grid <- cty %in% GULF_COUNTIES
    a <- if (in.grid) fn.assign_port(cty, CL, pop.center) else NULL
    if (in.grid && is.null(a)) {
      warning("No coastal land cell found in ", cty)
      in.grid <- FALSE
    }
    if (is.null(a)) a <- PORT_BLANK
    hd <- data.frame(County = cty, rank = which(bc$County == cty),
                     value = maj[[value.col]][maj$County == cty],
                     pct = maj$pct[maj$County == cty],
                     cum = maj$cum[maj$County == cty],
                     in_grid = in.grid, stringsAsFactors = FALSE)
    names(hd)[names(hd) == "value"] <- value.col
    res <- rbind(res, cbind(hd, a))
  }
  res[order(res$rank), ]
}


#' Write the binary port grid.
#'
#' 1 = port (a coastal land cell), 0 = water, NA/-9999 = all other land.
#'
#' @param depth SpatRaster template.
#' @param port.cells Cell numbers of the assigned ports.
#' @param file Output .asc path.
#' @return SpatRaster, invisibly.
fn.write_port_raster <- function(depth, port.cells, file) {
  dv <- terra::values(depth, mat = FALSE)
  # Every port must land on a NoData (land) cell of the depth grid.
  stopifnot(all(is.na(dv[port.cells])))

  out <- terra::rast(depth)
  terra::values(out) <- ifelse(is.na(dv), NA, 0)
  out[port.cells] <- 1

  if (!is.null(file)) {
    terra::writeRaster(out, file, overwrite = TRUE, NAflag = -9999,
                       datatype = "INT2S", filetype = "AAIGrid")
    unlink(paste0(file, ".aux.xml"))
  }
  invisible(out)
}


# ---- fleet definitions ------------------------------------------------------

#' Default fleet specifications.
#'
#' Each entry gives the data source, the row filter, the majority threshold and
#' the value column. Thresholds are carried over from the legacy scripts, where
#' they were tuned per fleet: broad fisheries concentrate in fewer counties, so
#' a lower threshold still captures the fleet.
#'
#' @return Named list of fleet specs.
fn.default_port_fleets <- function() {

  shrimp.spp <- c("SHRIMP, PINK", "SHRIMP, BROWN", "SHRIMP, WHITE",
                  "SHRIMP, ROCK", "SHRIMP, ROYAL RED", "SHRIMP, OTHER")
  reef.prefix <- c("GROUPER,", "SNAPPER,")
  reef.extra  <- c("GRUNTS", "AMBERJACKS")

  list(
    ports = list(
      source = "fwc", cum.thresh = 80, value.col = "Pounds",
      filter = function(d) d,
      label  = "All commercial landings, 1985-2025"),

    ports_grouper = list(
      source = "fwc", cum.thresh = 80, value.col = "Pounds",
      filter = function(d) d[d$Species %in% c("GROUPER, GAG", "GROUPER, RED"), ],
      label  = "Commercial gag + red grouper"),

    ports_gag = list(
      source = "fwc", cum.thresh = 90, value.col = "Pounds",
      filter = function(d) d[d$Species == "GROUPER, GAG", ],
      label  = "Commercial gag grouper"),

    ports_red = list(
      source = "fwc", cum.thresh = 90, value.col = "Pounds",
      filter = function(d) d[d$Species == "GROUPER, RED", ],
      label  = "Commercial red grouper"),

    ports_reef = list(
      source = "fwc", cum.thresh = 95, value.col = "Pounds",
      filter = function(d) {
        keep <- Reduce(`|`, lapply(reef.prefix, function(p) startsWith(d$Species, p))) |
                d$Species %in% reef.extra
        d[keep, ]
      },
      label  = "Commercial reef fish (groupers, snappers, grunts, amberjacks)"),

    ports_shrimp = list(
      source = "fwc", cum.thresh = 90, value.col = "Pounds",
      filter = function(d) d[d$Species %in% shrimp.spp, ],
      label  = "Commercial food shrimp (bait shrimp excluded)"),

    ports_baitshrimp = list(
      source = "fwc", cum.thresh = 90, value.col = "Pounds",
      filter = function(d) d[d$Species == "SHRIMP, BAIT", ],
      label  = "Commercial bait shrimp"),

    ports_headboat = list(
      source = "headboat", cum.thresh = 90, value.col = "Boats",
      pop.center = "headboat",
      label  = "Headboat / party boat fleet"),

    ports_rec_gag_charter = list(
      source = "mrip", cum.thresh = 90, value.col = "Dtrips",
      filter = function(m) m[m$common == "GAG" & m$mode_fx %in% 5, ],
      label  = "Gag, charter boat (MRIP mode_fx 5)"),

    ports_rec_gag_private = list(
      source = "mrip", cum.thresh = 90, value.col = "Dtrips",
      filter = function(m) m[m$common == "GAG" & m$mode_fx %in% 7, ],
      label  = "Gag, private/rental boat (MRIP mode_fx 7)"),

    ports_rec_red_rec = list(
      source = "mrip", cum.thresh = 90, value.col = "Dtrips",
      filter = function(m) m[m$common == "RED GROUPER" & m$mode_fx %in% c(5, 7), ],
      label  = "Red grouper, all recreational modes")
  )
}


#' Anchor table for a fleet: POP_CENTER, with headboat marinas swapped in.
fn.fleet_pop_center <- function(which = c("default", "headboat")) {
  which <- match.arg(which)
  pc <- POP_CENTER
  if (which == "headboat") {
    for (i in seq_len(nrow(HEADBOAT_HUB))) {
      j <- which(pc$County == HEADBOAT_HUB$County[i])
      if (length(j))
        pc[j, c("Town", "lon", "lat")] <- HEADBOAT_HUB[i, c("Town", "lon", "lat")]
    }
  }
  pc
}


# ---- entry point ------------------------------------------------------------

#' Build port basemaps for every fleet.
#'
#' @param depth SpatRaster of positive depths (the template grid).
#' @param dir.data Directory holding the landings and MRIP CSVs.
#' @param dir.maps Output directory.
#' @param fleets Fleet specs; see fn.default_port_fleets().
#' @param file.fwc,file.mrip Input filenames within `dir.data`.
#' @param write.ascii Write Ecospace ASCII grids and assignment CSVs.
#' @param verbose Report the majority set and port counties per fleet.
#' @return list(ports = SpatRaster one layer per fleet, assignments = named list
#'   of data.frames, coastal = the shared coastal-cell index).
fn.make_port_maps <- function(depth, dir.data, dir.maps,
                              fleets      = fn.default_port_fleets(),
                              file.fwc    = "ReportCreatorResults-County.csv",
                              file.mrip   = "MRIP WFS gag and red grouper dtrips by county.csv",
                              write.ascii = TRUE,
                              verbose     = TRUE) {

  if (inherits(depth, "Raster")) depth <- terra::rast(depth)
  stopifnot(inherits(depth, "SpatRaster"))
  if (is.na(terra::crs(depth)) || !nzchar(terra::crs(depth)))
    terra::crs(depth) <- "EPSG:4326"

  res.min <- round(terra::res(depth)[1] * 60, 0)
  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)

  # Shared spatial work, done once rather than per fleet.
  if (verbose) message("Building coastal candidate cells...")
  flc <- fn.florida_counties()
  CL  <- fn.build_coastal(depth, flc, verbose = verbose)

  needs <- vapply(fleets, function(f) f$source, character(1))
  fwc  <- if ("fwc"  %in% needs) fn.read_fwc_landings(file.path(dir.data, file.fwc))  else NULL
  mrip <- if ("mrip" %in% needs) fn.read_mrip_dtrips(file.path(dir.data, file.mrip))  else NULL

  lyrs <- list(); assigns <- list()

  for (nm in names(fleets)) {
    spec <- fleets[[nm]]
    if (verbose) message("---- ", nm, " ----")

    # Collapse the fleet's data to one value per county.
    bc <- switch(spec$source,
      fwc = {
        s <- spec$filter(fwc)
        stats::aggregate(stats::as.formula("Pounds ~ County"), s, sum)
      },
      mrip = {
        s <- spec$filter(mrip)
        a <- stats::aggregate(stats::as.formula("dtrip ~ County"), s, sum)
        names(a)[names(a) == "dtrip"] <- "Dtrips"
        a
      },
      headboat = HEADBOATS[HEADBOATS$Boats > 0, ],
      stop("Unknown source '", spec$source, "' for fleet ", nm))

    if (nrow(bc) == 0) {
      warning("Fleet '", nm, "' has no records after filtering; skipping.")
      next
    }

    bc  <- fn.county_majority(bc, spec$value.col, spec$cum.thresh)
    maj <- bc[bc$majority, ]
    if (verbose)
      message(sprintf("  majority (>= %g%%): top %d counties = %.1f%%",
                      spec$cum.thresh, attr(bc, "n_major"),
                      bc$cum[attr(bc, "n_major")]))

    pc  <- fn.fleet_pop_center(if (isTRUE(spec$pop.center == "headboat"))
                                 "headboat" else "default")
    res <- fn.port_table(maj, bc, CL, spec$value.col, pc)
    port <- res[res$in_grid, ]

    if (verbose) {
      message(sprintf("  in-grid port counties (%d): %s",
                      nrow(port), paste(port$County, collapse = ", ")))
      excl <- res$County[!res$in_grid]
      if (length(excl))
        message("  excluded (Atlantic / no in-grid coast): ",
                paste(excl, collapse = ", "))
      dup <- port$cell[duplicated(port$cell)]
      if (length(dup))
        message("  NOTE shared cells: ", paste(unique(dup), collapse = ", "))
    }

    f.asc <- if (write.ascii)
      file.path(dir.maps, sprintf("%s_%dmin.asc", nm, res.min)) else NULL
    r <- fn.write_port_raster(depth, port$cell, f.asc)
    names(r) <- nm
    lyrs[[nm]] <- r

    res$fleet <- nm
    assigns[[nm]] <- res
    if (write.ascii)
      utils::write.csv(res, file.path(dir.maps, sprintf("%s_assignment.csv", nm)),
                       row.names = FALSE)
  }

  if (length(lyrs) == 0) stop("No port layers were built.")
  ports <- terra::rast(lyrs)
  names(ports) <- names(lyrs)

  if (write.ascii) {
    all.assign <- do.call(rbind, lapply(assigns, function(a)
      a[, c("fleet", "County", "rank", "pct", "cum", "in_grid", "Town",
            "town_km", "nwater", "row", "col", "cell", "lon", "lat")]))
    utils::write.csv(all.assign,
                     file.path(dir.maps, sprintf("ports_all_assignments_%dmin.csv", res.min)),
                     row.names = FALSE)
    if (verbose) message("Ecospace ascii files written to\n  ", dir.maps)
  }

  invisible(list(ports = ports, assignments = assigns, coastal = CL))
}


# ---- plotting ---------------------------------------------------------------

#' Validation map for one fleet: port cells on the depth grid.
#'
#' @param depth SpatRaster template.
#' @param port In-grid rows of the assignment table.
#' @param CL Output of fn.build_coastal().
#' @param pop.center Anchor table used for this fleet.
#' @param file Output PNG path.
#' @param main Plot title.
#' @return The path written, invisibly.
fn.plot_port_validation <- function(depth, port, CL, pop.center = POP_CENTER,
                                    file, main = "") {

  brks <- c(seq(0, 200, 20), seq(300, 1000, 100), 3000)
  dmax <- max(terra::values(depth, mat = FALSE), na.rm = TRUE)
  brks <- c(brks[brks < dmax], ceiling(dmax))
  cols <- colorRampPalette(c("#cfe8ff", "#4a90d9", "#0b3d91"))(length(brks) - 1)

  png(file, width = 8, height = 8.5, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)

  terra::plot(depth, col = cols, breaks = brks, colNA = "grey80",
              mar = c(4, 4, 3, 6), plg = list(cex = 0.7), main = main)
  maps::map("county", "florida", add = TRUE, col = "grey40", lwd = 0.4)

  points(terra::xyFromCell(depth, CL$coast_cell), pch = 15, col = "grey60", cex = 0.45)

  if (nrow(port) > 0) {
    pcxy <- terra::xyFromCell(depth, port$cell)
    pt   <- pop.center[match(port$County, pop.center$County), ]
    segments(pt$lon, pt$lat, pcxy[, 1], pcxy[, 2], col = "red", lwd = 0.8)
    points(pt$lon, pt$lat, pch = 21, bg = "red", col = "black", cex = 0.9)
    points(pcxy, pch = 22, bg = "yellow", col = "black", cex = 1.7, lwd = 1.2)
    text(pcxy[, 1], pcxy[, 2], port$County, pos = 3, cex = 0.6, font = 2)
  }
  legend("bottomleft", bty = "n", cex = 0.6, pt.cex = c(1.4, 0.9, 0.6),
         pch = c(22, 21, 15), pt.bg = c("yellow", "red", NA),
         col = c("black", "black", "grey60"),
         legend = c("port cell (coastal land)", "anchor point",
                    "coastal land candidates"))
  box()
  invisible(file)
}


#' Write a validation map for every fleet, plus a panel of all port layers.
#'
#' @param x list from fn.make_port_maps().
#' @param depth SpatRaster template.
#' @param dir.maps Output directory.
#' @param fleets Fleet specs, used for the titles.
#' @return Character vector of files written, invisibly.
fn.plot_port_maps <- function(x, depth, dir.maps,
                              fleets = fn.default_port_fleets()) {

  if (!dir.exists(dir.maps)) dir.create(dir.maps, recursive = TRUE)
  res.min <- round(terra::res(depth)[1] * 60, 0)
  out <- character(0)

  for (nm in names(x$assignments)) {
    res  <- x$assignments[[nm]]
    port <- res[res$in_grid, ]
    spec <- fleets[[nm]]
    pc <- fn.fleet_pop_center(if (!is.null(spec$pop.center) &&
                                  isTRUE(spec$pop.center == "headboat"))
                                "headboat" else "default")
    f <- file.path(dir.maps, sprintf("%s_validation_%dmin.png", nm, res.min))
    fn.plot_port_validation(depth, port, x$coastal, pc, f,
                            main = sprintf("%s\n%d port cells",
                                           if (!is.null(spec$label)) spec$label else nm,
                                           nrow(port)))
    out <- c(out, f)
  }

  # Panel of every fleet's port cells.
  p <- x$ports
  n <- terra::nlyr(p)
  nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)
  fp <- file.path(dir.maps, sprintf("ports %dmin.png", res.min))
  png(fp, height = 3.2 * nr, width = 3.6 * nc, units = "in", res = 300)
  par(mfrow = c(nr, nc))
  for (i in seq_len(n)) {
    v <- terra::values(p[[i]], mat = FALSE)
    n.port <- sum(v == 1, na.rm = TRUE)
    terra::plot(p[[i]], col = c("gray95", "firebrick"), breaks = c(-0.5, 0.5, 1.5),
                colNA = "lightgray", legend = FALSE, mar = c(2, 2, 3.5, 1),
                main = sprintf("%s\n%d port cells", names(p)[i], n.port))
    maps::map("state", add = TRUE, col = "gray40")
  }
  dev.off()
  out <- c(out, fp)

  message("Figures written to\n  ", dir.maps)
  invisible(out)
}
