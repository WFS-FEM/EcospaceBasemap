#Control script for the WFS Ecospace basemaps.
#
#Every path below is repo-relative, so a fresh clone runs without editing this
#file. To read data from, or write outputs to, somewhere else, put a
#config.local.R in the repo root - it is gitignored, and config.local.example.R
#shows what goes in it.
#
#Written to be stepped through interactively: each section leaves its result in
#the workspace (depth, seagrass, seabed, gfisher, ar, basemap, ma, ports,
#regions) so you can inspect or re-plot without re-running the expensive parts.

rm(list=ls());graphics.off();gc()
if(interactive() && .Platform$OS.type=='windows'){
  try(rm(.SavedPlots,envir=.GlobalEnv),silent=TRUE)  #RGui plot history
  try(dev.new(record=TRUE),silent=TRUE)
}

#everything downstream is resolved from the working directory, so check it first
if(!file.exists('EcospaceBasemap.Rproj'))
  stop('Set the working directory to the repo root - open EcospaceBasemap.Rproj, ',
       'or setwd() there. Currently: ',getwd())

#load (source) all the functions in the R folder
#sourced BEFORE the library() calls below: if a function file ever attaches
#raster or sp, doing it after library(terra) would mask terra's generics
#(rasterize, crs, resample, extract) for the whole session.
invisible(sapply(list.files(file.path(getwd(),"R"),full.names=T),source))

library('marmap')
library('terra')
library('sf')
library('colorRamps')


#SETUP---------------------
res <- 5 #minutes
bbox <- c(-87.5,-81,25,30.5) #decimal degrees
excl.depth <- 500 #meters

#paths - repo-relative defaults, overridable in config.local.R
dir.data     <- file.path(getwd(),'data')                     #source data
dir.basemaps <- file.path(getwd(),'output',paste0(res,'min')) #generated grids
file.gdb     <- NULL  #NULL discovers a single .gdb inside dir.data

if(file.exists('config.local.R')){
  source('config.local.R')
  message('Applied local path overrides from config.local.R')
}

dir.depth <- file.path(dir.basemaps,'depth')
dir.habitats <- file.path(dir.basemaps,"habitat")
invisible(lapply(c(dir.basemaps,dir.depth,dir.habitats),dir.create,
                 recursive=TRUE,showWarnings=FALSE))

#INPUT DATA---------------------
#Most of data/ is gitignored, so a clone starts with 3 of the 9 inputs present.
#fn.pull_all() fetches everything that downloads itself and skips what is already
#there; fn.check_inputs() then reports what is present, stops if a REQUIRED input
#is missing, and prints where to obtain it. R/data_setup_functions.R holds the
#manifest both of them read from.
fn.pull_all(dir.data)
have <- fn.check_inputs(dir.data, file.gdb=file.gdb)

#1. depth-----------
##pull data----
depth = getNOAA.bathy(lon1=bbox[1],lon2=bbox[2],lat1=bbox[3],lat2=bbox[4],resolution=res)
depth = marmap::as.raster(depth)
depth[depth>0] = NA
depth = depth*-1
depth[depth==0] = min(depth[depth>0])
xy <- xyFromCell(depth, 1:2)
dist <- round(as.numeric(terra::distance(xy, lonlat = TRUE, unit = "km")), 1)
depth <- rast(depth)

##exclusion layer----
excl <- depth
excl[depth>excl.depth] = 1
excl[depth<=excl.depth] = NA

##output ascii----
file.depth <- file.path(dir.depth,paste0('depth_',res,'min.asc'))
file.excl <- file.path(dir.depth,paste0('excl_',res,'min.asc'))
terra::writeRaster(depth,file.depth, overwrite=T)
terra::writeRaster(excl,file.excl, overwrite=T)

##plots----
colv    = c('lightblue','blue','darkblue')
funpal  = colorRampPalette(colv,bias=2)
brks    = c(seq(0,100,20),200,300,400,500,1000,terra::minmax(depth)["max", ])
nbcols  = length(brks)-1
color.dep   = funpal(nbcols)

png(file.path(dir.depth,paste0('depth_',res,'min.png')),width = 6.5, height=9, units='in',res=300)
plot(depth,colNA='black',main=paste0(round(res(depth)*60)[1],' min / ',dist,' km / ',dim(depth)[1],'x',dim(depth)[2]),col=color.dep,breaks=brks)
dev.off()  

#2. habitats--------
##2.1 seagrass--------------------------------------------------------------------------------------
###output----
dir.create(file.path(dir.habitats,'seagrass'), recursive=T, showWarnings=F)
seagrass.prop.gulf <- fn.make_seagrass_ascii(dir.seagrass= file.path(dir.data,'seagrass',"GulfwideSAV"),
                                             dir.ascii = file.path(dir.habitats,'seagrass'), depth = depth)
plot(seagrass.prop.gulf)

if(have[['seagrass_fwc']]){
  seagrass.prop.fwc <- fn.make_seagrass_ascii(dir.seagrass= file.path(dir.data,'seagrass',"Seagrass_Statewide"),
                                              dir.ascii = file.path(dir.habitats,'seagrass'), depth = depth)
  plot(seagrass.prop.fwc)
}

###combine the two sources----
#The two layers map the SAME beds, not different ones (186 of 247 Gulfwide cells
#and 186 of 212 FWC cells overlap), so they cannot be added - summing would push
#139 cells past 1. This single combined map is what feeds the sum-to-1 step.
#
#max() is the LOWER bound of the union: for a cell holding a and b the truth is
#a + b - overlap, and the rasters don't record the overlap. It is exact here
#because FWC beds are nested inside Gulfwide beds, so overlap = min(a,b) and
#a + b - min(a,b) collapses to max(a,b) - verified against the polygon route at
#5 min, which gave an identical map. Re-check with fn.combine_seagrass() if you
#build at a different resolution, since nesting is a property of the polygons
#and the cell size sets how much mixing happens within a cell.
#
#Vintage caveat: GulfwideSAV carries hab_88/hab_92 fields and reports ~2x the
#coverage of the newer FWC layer, so the union is closer to historical maximum
#extent than to current extent.
#
#Without the FWC layer there is nothing to combine: GulfwideSAV is the wider of
#the two anyway (roughly twice the coverage), so the map is the union minus
#whatever FWC maps and Gulfwide does not. Section 2.5 is unaffected in shape.
if(have[['seagrass_fwc']]){
  seagrass <- fn.combine_seagrass_rasters(c(seagrass.prop.gulf, seagrass.prop.fwc),
                                          method    = 'max',
                                          depth     = depth,
                                          dir.ascii = file.path(dir.habitats,'seagrass'),
                                          label     = 'combined')
} else {
  warning('FWC Seagrass_Statewide is missing - section 2.1 uses GulfwideSAV alone. ',
          'Run fn.pull_all(dir.data) to fetch it.', call.=FALSE)
  seagrass <- seagrass.prop.gulf
}
plot(seagrass)

#the polygon route, for validation - slow (~5 min: reading and repairing ~90k
#features dominates, not the rasterizing):
#seagrass <- fn.combine_seagrass(dirs.seagrass = file.path(dir.data,'seagrass',
#                                                          c("GulfwideSAV","Seagrass_Statewide")),
#                                dir.ascii = file.path(dir.habitats,'seagrass'),
#                                depth = depth)

##2.2 dbSeabed----
dir.dbseabed <- file.path(dir.data,'dbseabed')

###make rasters----
seabed <- fn.rasterize_dbseabed(dir.dbseabed = file.path(dir.data,'dbseabed'), dir.out=file.path(dir.habitats,'dbseabed'),
                                depth, resample.method = 'near')

###plots----
#grey patches are NA gaps in dbSEABED; section 2.5 fills them with a focal mean
fn.plot_dbseabed(seabed, dir.maps = file.path(dir.habitats,'dbseabed'))

##2.3 GFISHER----
dir.create(file.path(dir.habitats,'gfisher'), recursive = TRUE, showWarnings = FALSE)
#set explicitly with file.gdb in config.local.R, or discovered by extension
file.gfishergdb = if(!is.null(file.gdb)) file.gdb else
  list.files(dir.data,pattern="\\.gdb$",full.names = T)[1]

###make rasters----
#natural classes are extrapolated beyond the side-scan footprint out to the
#exclusion depth; artificial classes are left at 0 where unmapped.
#fill.method='idw' is the legacy distance-only fill (idp=4, nmax=8), 
#'gam' (predicts from depth + dbSEABED substrate, pass covariates=seabed) and
#'strata' (depth-bin ratio estimator) are the alternatives.
gfisher <- fn.make_GFISHER_habitat_maps(depth         = depth,
                                        file.gdb      = file.gfishergdb,
                                        dir.maps      = file.path(dir.habitats,'gfisher'),
                                        max.depth.hab = 200,
                                        fill.method   = 'idw',
                                        anchor.zero = 'both')

###plots----
fn.plot_GFISHER_habitats(gfisher, dir.maps=file.path(dir.habitats,'gfisher'), 
                         col = colorRampPalette(colorRamps::matlab.like(100), bias = 3)(100))


##2.4 artificial reefs----
#Structure records (dataS2) carry footprint area but no relief; the FWC
#deployment table (reeflocations) carries relief but is keyed only by free-text
#description, so the two are joined by fuzzy string matching. Slow (~2 min at
#5 min resolution) - the fuzzy joins dominate, not the rasterizing.
#Without reeflocations.csv there is no relief, so the run falls back to
#unweighted footprint area and writes the Low/Medium/High layers as zeros.
dir.ar <- file.path(dir.habitats,'artificial_reefs')
dir.create(dir.ar, recursive=T, showWarnings=F)

ar <- fn.make_AR_maps(depth     = depth,
                      file.ar   = file.path(dir.data,'artificial_reefs',
                                            'dataS2_artificial_reef_structures_REDACTED.csv'),
                      file.reef = if(have[['artificial_reefs_fwc']])
                                    file.path(dir.data,'artificial_reefs','reeflocations.csv'),
                      dir.maps  = dir.ar)

###plots----
fn.plot_AR_maps(ar, dir.maps = dir.ar)


##2.5 combine to sum 1----
#Ten layers that sum to exactly 1 in every water cell, which is what Ecospace
#needs of a habitat basemap:
#  AL AM AH   artificial reef (GFISHER + the AR database, matched by name)
#  NL NM NH   natural reef    (GFISHER + dbSEABED rock, dissolved by stratum)
#  seagrass
#  Sand Mud Gravel
#Rock is split across the three natural classes using the relief composition
#observed in mapped cells within each depth x 1-degree-latitude stratum, falling
#back to depth bin then global where a stratum is thin. RCK_raw is also written,
#outside the sum, so the undissolved layer can still be inspected.
dir.sum1 <- file.path(dir.habitats,'sum1')
dir.create(dir.sum1, recursive=T)

basemap <- fn.combine_habitats_sum1(depth     = depth,
                                    gfisher   = gfisher$habitat,
                                    ar        = ar,
                                    seabed    = seabed,
                                    seagrass  = seagrass,
                                    microgrid = gfisher$microgrid,
                                    dir.maps  = dir.sum1)

###plots----
fn.plot_habitat_basemap(basemap, dir.maps = dir.sum1)


#3. management areas------------
#One grid per management area, aligned to the depth template:
#  1      cells the area covers
#  0      other cells inside the model domain
#  -9999  land and excluded deep water
#touches=TRUE so small or thin areas are not lost at coarse resolution.
#madswan_steamboat_edges holds three features in a LABEL field and is split into
#two grids - Madison/Swanson + Steamboat Lumps share seasonal management, the
#Edges differ. See fn.default_ma_splits() to add or change split rules.
dir.ma <- file.path(dir.basemaps,'management_areas')
dir.create(dir.ma, recursive=T)

ma <- fn.make_management_area_maps(depth    = depth,
                                   dir.zips = file.path(dir.data,'management_areas'),
                                   dir.maps = dir.ma)

###plots----
fn.plot_management_areas(ma, dir.maps = dir.ma)


#4. ports-----------
#One binary grid per fleet: 1 = port, 0 = water, -9999 = all other land.
#A port is a COASTAL LAND cell of the depth grid (land = NoData in depth,
#touching >=1 water cell), inside the county, nearest that county's anchor
#point. Ports go to the smallest set of counties reaching the fleet's
#cum.thresh share of landings / trips / vessels.
#Gulf vs Atlantic is set by the explicit GULF_COUNTIES list, not by distance -
#the grid's east edge clips Atlantic water near Jacksonville and Cape Canaveral.
#11 fleets; edit fn.default_port_fleets() to change filters or thresholds, and
#POP_CENTER / HEADBOAT_HUB in R/port_functions.R to move a port.
dir.ports <- file.path(dir.basemaps,'ports')
dir.create(dir.ports, recursive=T)

ports <- fn.make_port_maps(depth    = depth,
                           dir.data = file.path(dir.data,'ports'),
                           dir.maps = dir.ports)

###plots----
#one validation map per fleet, plus a panel of all 11 port layers
fn.plot_port_maps(ports, depth = depth, dir.maps = dir.ports)


#5. regions---------
#One categorical grid coding which survey region each model cell belongs to:
#  -9999  land
#      0  water, unsampled
#   1-9   age-0 survey regions (bays and estuaries)
#  10-16  GFISHER survey coverage: 10 FWRI, 11 PASC, 12 PC,
#         13 F+P, 14 F+PC, 15 P+PC, 16 F+P+PC
#Age-0 regions override GFISHER codes where they overlap.
#
#The digitizing step is NOT reproduced here. It extracted the age-0 polygons
#from a georeferenced PNG, and R's PNG decoder renders the anti-aliased borders
#about a pixel thinner than Python's, which can resolve only 8 of the 9 regions.
#age0_survey_regions.shp from the Python pipeline is the authoritative input.
dir.regions <- file.path(dir.basemaps,'regions')
dir.create(dir.regions, recursive=T)
dir.regdata <- file.path(dir.data,'regions')

###age-0 regions----
#regenerates the UN-edited grid; the combine step below uses the hand-edited one
age0 <- fn.make_age0_region_grids(depth,
                                  file.shp = file.path(dir.regdata,'age0_survey_regions.shp'),
                                  dir.maps = dir.regions)

###GFISHER survey coverage----
#kernel density on the sample points, keep the densest 95%, concave hull.
#O(n^2) KDE to match scipy's gaussian_kde - takes ~30 s at 5 min.
gfisher.reg <- fn.make_gfisher_survey_regions(depth,
                                              file.csv = file.path(dir.regdata,'env3LABS_93to24.csv'),
                                              dir.maps = dir.regions)

###combine----
#age0_survey_regions_5min_mod.asc is HAND-EDITED and is an input, not a
#derived product - pass the regenerated `age0` instead to skip those edits.
regions <- fn.combine_regions(age0    = file.path(dir.regdata,'age0_survey_regions_5min_mod.asc'),
                              gfisher = gfisher.reg,
                              dir.maps = dir.regions)

###attributes + plot----
region.attr <- fn.region_attributes(regions, dir.maps = dir.regions)
fn.plot_regions(regions, dir.maps = dir.regions)




