#Local path overrides for make_WFS_basemaps.R.
#
#Copy this file to config.local.R in the repo root and edit it. config.local.R is
#gitignored, so machine-specific paths never reach the repository - this is the
#only file you should need to touch to run the pipeline somewhere else.
#
#The driver sources it AFTER setting its repo-relative defaults and BEFORE using
#any of them, so uncomment only the lines you actually want to change. `res` is
#already defined at that point, which is why it can be used below.


#Where the source data lives.
#Default: file.path(getwd(),'data')
#Point this at a shared drive to avoid keeping a second copy of the 600 MB of
#inputs, or to let several clones share one data store.
#dir.data <- 'D:/shared/WFS/EcospaceBasemap-data'


#Where the generated grids are written.
#Default: file.path(getwd(),'output',paste0(res,'min'))
#dir.basemaps <- file.path('C:/Users/you/OneDrive - University of Florida',
#                          'WFS Fisheries Ecosystem Modeling/WFS EwE/Ecospace/basemaps',
#                          paste0(res,'min'))


#The GFISHER geodatabase, if you do not want a 274 MB copy inside data/.
#Default: NULL, which discovers a single .gdb inside dir.data.
#file.gdb <- file.path('C:/Users/you/University of Florida',
#                      'Chagaris, David - WFS Fisheries Ecosystem Modeling',
#                      'data/GFISHER/April2026/GFISHER_EAST_Universe_2026.gdb')
