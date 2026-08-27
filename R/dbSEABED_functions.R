library('terra')


#' Download and Extract dbSEABED Grain-Size Layers for the Northern Gulf of Mexico
#'
#' This function downloads the dbSEABED ZIP file sets for
#' percent gravel, sand, mud, and rock for the Northern Gulf of Mexico
#' from the CSDMS dbSEABED data repository and extracts them to a local directory.
#'
#' @details
#' The function retrieves four ZIP archives corresponding to:
#' \itemize{
#'   \item \code{Gmf_GVL} – percent gravel
#'   \item \code{Gmf_SND} – percent sand
#'   \item \code{Gmf_MUD} – percent mud
#'   \item \code{Gmf_RCK} – percent rock
#' }
#'
#' Each ZIP file is downloaded to the specified output directory,
#' extracted into its own subfolder (named after the ZIP filename),
#' and the ZIP file is deleted after extraction.
#'
#' Source data and description:
#' \url{https://csdms.colorado.edu/wiki/DBSEABED#Data_for_Modellers}
#'
#' @param dir.out Character string giving the directory where all ZIP files
#'   will be downloaded and extracted. Subdirectories will be created automatically
#'   based on the ZIP file names.
#'
#' @return
#' Invisibly returns \code{NULL}. Files are written to \code{dir.out} as a side effect.
#'
#' @examples
#' \dontrun{
#' fn.pull_dbseabed(
#'   dir.out = "C:/path/to/dbSEABED"
#' )
#' }
#'
#' @export
fn.pull_dbseabed <- function(dir.out){
  
  #dir.out = "C:\\Users\\dchagaris\\Github\\WFS-FEM\\EnvironmentalDrivers2EwE\\data\\dbSEABED"
  
  message('dbSEABED Zipfile sets for % gravel, sand, mud, and rock in Northern Gulf of Mexico will be downloaded from\nhttps://csdms.colorado.edu/wiki/DBSEABED#Data_for_Modellers') 
  
  url.set <- c("https://csdms.colorado.edu/csdms_wiki/images/Gmf_GVL.zip",
                "https://csdms.colorado.edu/csdms_wiki/images/Gmf_SND.zip",
                "https://csdms.colorado.edu/csdms_wiki/images/Gmf_MUD.zip",
                "https://csdms.colorado.edu/csdms_wiki/images/Gmf_RCK.zip")

  
  for(i in 1:length(url.set)){
    url = url.set[i]
    file.out = file.path(dir.out,basename(url))
    download.file(url, destfile = file.out, mode='wb')
    unzip(zipfile=file.out, exdir=file.path(dir.out, gsub(".zip","",basename(url))))
    unlink(file.out)
  }
  
  message('All layers downloaded and extracted to \n',dir.out)
} #eof


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#' Process and Convert dbSEABED Grain-Size Layers to Aligned ASCII Grids
#'
#' This function reads raw dbSEABED percent-composition layers
#' (gravel, sand, mud, rock) stored as ESRI ASCII files, harmonizes each layer
#' to match a target template grid, rescales the data to proportions,
#' renormalizes layers so they sum to 1, and writes each output map to an
#' ASCII grid suitable for modeling applications (e.g., EwE/Ecospace).
#'
#' @details
#' For each subdirectory inside \code{dir.dbseabed}, the function:
#' \enumerate{
#'   \item Loads the source ASCII raster(s).
#'   \item Ensures the raster has the same CRS as \code{depth}; reprojects if necessary.
#'   \item Crops to the template grid extent.
#'   \item Zeroes negative values and rescales percentages to proportions.
#'   \item Aggregates to approximately match the resolution of \code{depth}.
#'   \item Resamples onto the \code{depth} grid using the specified method.
#'   \item Collects all layers into a single multi-layer object.
#' }
#'
#' After processing all layers, the function normalizes the resulting rasters
#' so that the four dbSEABED fractions sum to 1 at every grid cell.  
#' Each normalized layer is then exported as an ESRI ASCII grid
#' (\code{.asc}) into \code{dir.ascii}.
#'
#' A PDF file showing the processing stages (original, cropped, aggregated,
#' and final aligned raster) is created as a diagnostic output.
#'
#' @param dir.dbseabed Character string giving the path to the directory
#'   containing the extracted dbSEABED folders (e.g., \code{"Gmf_GVL"}, \code{"Gmf_SND"}, etc.).
#' @param dir.ascii Character string giving the directory where output ASCII files
#'   and diagnostic plots will be written. The directory is created if it does not exist.
#' @param depth A \code{SpatRaster} (or legacy \code{RasterLayer}) used as the
#'   template grid for final alignment. Its resolution, extent, and CRS determine
#'   the target output grid.
#' @param resample.method Character string giving the interpolation method passed
#'   to \code{terra::resample()}. Typical options include \code{"near"} (default)
#'   and \code{"bilinear"}.
#'
#' @return
#' Invisibly returns \code{NULL}.  
#' The function writes:
#' \itemize{
#'   \item a set of normalized ESRI ASCII files to \code{dir.ascii}, and
#'   \item a PDF file (\code{dbsEABED_plots.pdf}) showing intermediate processing steps.
#' }
#'
#' @examples
#' \dontrun{
#' fn.make_dbseabed_ascii(
#'   dir.dbseabed = "C:/path/to/dbSEABED",
#'   dir.ascii    = "C:/path/to/ascii_outputs",
#'   depth        = depth.15min,
#'   resample.method = "near"
#' )
#' }
#'
#' @export
fn.make_dbseabed_ascii <- function(dir.dbseabed, dir.ascii, depth, resample.method='near'){
  
  # dir.dbseabed = "C:\\Users\\dchagaris\\Github\\WFS-FEM\\EnvironmentalDrivers2EwE\\data\\dbSEABED"
  # dir.ascii = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\WFS Fisheries Ecosystem Modeling\\WFS EwE\\Ecospace\\maps\\dbsEABED"
  # depth = terra::rast(depth.15min)
  
  dirs.dbseabed = list.dirs(dir.dbseabed, recursive=F)
  
  res <- res(depth)[1]*60
  
  
  # Ensure 'depth' is a SpatRaster (terra) not a RasterLayer (raster)  <-- added
  if (inherits(depth, "Raster")) {
    depth <- terra::rast(depth)
  }
  
  
  seabed.stack <- terra::rast()
  pdf(file=file.path(dir.ascii,'dbsEABED_plots.pdf'), onefile = T)
  for(i in 1:length(dirs.dbseabed)){
    #i=1
    file.asc.i = list.files(dirs.dbseabed[i],pattern=".asc$", full.names=T)
    rast.i = terra::rast(file.asc.i)
    #crs(rast.i)==crs(depth)
    
    
    ## ---- CRS harmonization -------------------------------------------------
    # Ensure CRS matches 'depth'; project if needed (minimal change)
    if (!terra::same.crs(rast.i, depth)) {  # <-- added
      rast.i <- terra::project(rast.i, depth, method = resample.method)  # <-- added
    }
    
    
    crop.i = terra::crop(rast.i,depth)
    crop.i[crop.i < 0] <- 0
    crop.i <- crop.i / 100

    approx_fact <- round(terra::res(depth) / terra::res(crop.i))
    approx_fact[approx_fact < 1] <- 1
    
    agg.i <- terra::aggregate(crop.i, fact = approx_fact, fun = mean, na.rm = TRUE)
    
    # 2) Then resample onto the target grid
    # Continuous:
    out.i <- terra::resample(agg.i, depth, method = resample.method)  # avoids extra smoothing
    out.i[is.na(depth)] <- NA

    
    seabed.stack <- c(seabed.stack, out.i)
    
    par(mfrow=c(2,2),oma=c(0,0,1,0))
    plot(rast.i)
    plot(crop.i)
    plot(agg.i)
    plot(out.i, colNA='gray')
    title(main=paste("dbSEABED:",strsplit(basename(file.asc.i),"_")[[1]][2]), outer=T)
  }
  dev.off()
  
  
  # Renormalize to sum 1
  sum_layers <- sum(seabed.stack, na.rm = TRUE)
  stack_norm <- (seabed.stack/sum_layers)
  #save
  for(i in 1:nlyr(stack_norm)){
    terra::writeRaster(stack_norm[[i]], filename=file.path(dir.ascii,paste0(gsub("val","prop",names(stack_norm)[i]),"_",res,"min.asc")),
                       overwrite=T, gdal=c('DECIMAL_PRECISION=4'), NAflag=-9999)
  }
  message('dbSEABED ascii files saved to\n',dir.ascii)

  #rm.xml <- list.files(dir.ascii, pattern=".xml$", full.names=T)
  #rm.prj <- list.files(dir.ascii, pattern=".prj$", full.names=T)

  
} #eof


fn.rasterize_dbseabed <- function(dir.dbseabed, depth, resample.method='near', dir.out){
  
  # dir.dbseabed = "C:\\Users\\dchagaris\\Github\\WFS-FEM\\EnvironmentalDrivers2EwE\\data\\dbSEABED"
  # dir.ascii = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\WFS Fisheries Ecosystem Modeling\\WFS EwE\\Ecospace\\maps\\dbsEABED"
  # depth = terra::rast(depth.15min)
  
  dirs.dbseabed = list.dirs(dir.dbseabed, recursive=F)
  dir.create(dir.out,recursive = T)
  res <- res(depth)[1]*60
  
  
  # Ensure 'depth' is a SpatRaster (terra) not a RasterLayer (raster)  <-- added
  if (inherits(depth, "Raster")) {
    depth <- terra::rast(depth)
  }
  
  
  seabed.stack <- terra::rast()
  pdf(file=file.path(dir.out,'dbseabed_plots.pdf'), onefile = T)
  for(i in 1:length(dirs.dbseabed)){
    #i=1
    file.asc.i = list.files(dirs.dbseabed[i],pattern=".asc$", full.names=T)
    rast.i = terra::rast(file.asc.i)
    #crs(rast.i)==crs(depth)
    
    ## ---- CRS harmonization -------------------------------------------------
    # Ensure CRS matches 'depth'; project if needed (minimal change)
    if (!terra::same.crs(rast.i, depth)) {  # <-- added
      rast.i <- terra::project(rast.i, depth, method = resample.method)  # <-- added
    }
    
    crop.i = terra::crop(rast.i,depth)
    crop.i
    crop.i[crop.i == -99] <- NA
    crop.i <- crop.i / 100
    approx_fact <- round(terra::res(depth) / terra::res(crop.i))
    approx_fact[approx_fact < 1] <- 1
    
    agg.i <- terra::aggregate(crop.i, fact = approx_fact, fun = mean, na.rm = TRUE)

    # 2) Then resample onto the target grid
    # Continuous:
    out.i <- terra::resample(agg.i, depth, method = resample.method)  # avoids extra smoothing
    out.i[is.na(depth)] <- NA
    seabed.stack <- c(seabed.stack, out.i)
    
    par(mfrow=c(2,2),oma=c(0,0,1,0))
    plot(rast.i, main='full raster')
    plot(crop.i, main='WFS cropped')
    plot(agg.i, main='Aggregated to approximate depth res')
    plot(out.i, main='Resampled to depth grid')
    title(main=paste("dbSEABED:",strsplit(basename(file.asc.i),"_")[[1]][2]), outer=T)
  }
  dev.off()
  writeRaster(seabed.stack,filename=paste0(dir.out,"/",names(seabed.stack),"_",res,"min.asc"))
return(seabed.stack)
} #eof


