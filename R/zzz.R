### use onLoad if I ever run Python through reticulate for the final analyses. Specify the packages I want to delay loading here, so users can load their own virtual environments for python.

### global reference to scipy (will be initialized in .onLoad)
# scipy <- NULL
#
# .onLoad <- function(libname, pkgname) {
#   # use superassignment to update global reference to scipy
#   scipy <<- reticulate::import("scipy", delay_load = TRUE)
# }
