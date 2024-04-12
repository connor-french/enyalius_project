# project models

# load packages
library(terra)
library(stringr)
library(purrr)
library(here)
library(ENMeval)
library(sf)



# define functions --------------------------------------------------------

## bio is a vector of bioclims to read in. It should only contain the bioclim name, e.g. "bio02"
## timeslice is which time period to read in for historical climate. Should be a single number
## when rename is TRUE, it renames the variables to agree with the variable names used in modeling
read_predictors <- function(timeslice, rename = TRUE) {

  all_paths <- list.files(here("analysis", "output", "historical_climate"), pattern = ".tif", full.names = TRUE)

  # paths that contain the bioclims for the selected time
  time_paths <- all_paths[str_detect(all_paths, paste0("_", timeslice, ".tif"))]

  # remove time slice from variable name
  if (timeslice == 0) {
    varnames_clean <- gsub(regex("_[0-9]+.tif"), "", time_paths) %>% basename()
  } else {
    varnames_clean <- gsub(regex("_.[0-9]+.tif"), "", time_paths) %>% basename()
  }


  # rename to align with the variables I used with the SDM models
  if (rename) {
    varnames_clean <- str_replace(varnames_clean, "bio0", "bio")
  }

  r_stack <- terra::rast(time_paths)

  names(r_stack) <- varnames_clean

  return(r_stack)

}

# I have to crop according to the bioclims because they have the coastline
crop_mask <- function(pred_rast, bioclim) {
  # crop by predictor
  c <- crop(pred_rast, bioclim)

  # mask
  m <- mask(c, bioclim)

  # add the predictor raster layer names back
  names(m) <- names(pred_rast)

  return(m)
}



# read in predictors ------------------------------------------------------

# read in predictors
preds <- lapply(rev(seq(-200, 20, 10)), read_predictors)

# make names of predictor timesteps reflect the year annotations
names(preds) <- paste0("y", 0:22, "k")

# read in cropped bioclims to crop and mask the predictors with
bioclims <- rast(here("analysis", "output", "cropped_predictors", "bioclims.tif"))

# mask predictor rasters by bioclims
preds_masked <- map(names(preds), \(x) crop_mask(pred_rast = preds[[x]], bioclim = bioclims$bio1))

# name according to predictor raster list names
names(preds_masked) <- names(preds)


# iheringii ---------------------------------------------------------------

# read in model
sdm_iheringii <- readRDS(here("analysis", "output", "sdm_models", "sdm_iheringii.rds"))

# predict to the full extent
projections_ihe <- map(preds_masked, \(x) dismo::predict(sdm_iheringii@models[["rm.2_fc.LQ"]], x))

# read in iheringii locs
locs_ihe <- st_read(here("analysis", "output", "thinned_localities", "iheringii_thinned.gpkg"))

# crop to iheringii locs with buffer
## larger buffer than used to train SDM to allow for dispersal
mcp_ihe <- st_convex_hull(st_union(locs_ihe)) %>%
  st_buffer(dist = units::set_units(2.0, degree)) %>%
  terra::vect()

projections_ihe_crop <- crop(rast(projections_ihe), mcp_ihe) %>%
  mask(mcp_ihe)

# plot stability of environment through time
plot(stdev(projections_ihe_crop, na.rm = TRUE))
points(vect(locs_ihe))

# write predictions to file
writeRaster(projections_ihe_crop, here("analysis", "output", "sdm_projections", "projections_ihe.tif"))


# catenatus ---------------------------------------------------------------

# read in model
sdm_catenatus <- readRDS(here("analysis", "output", "sdm_models", "sdm_catenatus.rds"))

# predict to the full extent
projections_cat <- map(preds_masked, \(x) dismo::predict(sdm_catenatus@models[["rm.0.5_fc.L"]], x))

# read in catenatus locs
locs_cat <- st_read(here("analysis", "output", "thinned_localities", "catenatus_thinned.gpkg"))

# crop to catenatus locs with buffer
## larger buffer than used to train SDM to allow for dispersal
mcp_cat <- st_convex_hull(st_union(locs_cat)) %>%
  st_buffer(dist = units::set_units(2.0, degree)) %>%
  terra::vect()

projections_cat_crop <- crop(rast(projections_cat), mcp_cat) %>%
  mask(mcp_cat)

# plot stability of environment through time
plot(stdev(projections_cat_crop, na.rm = TRUE))
points(vect(locs_cat))

# write predictions to file
writeRaster(projections_cat_crop, here("analysis", "output", "sdm_projections", "projections_cat.tif"))


# perditus ----------------------------------------------------------------
# read in model
sdm_perditus <- readRDS(here("analysis", "output", "sdm_models", "sdm_perditus.rds"))

# predict to the full extent
projections_per <- map(preds_masked, \(x) dismo::predict(sdm_perditus@models[["rm.2_fc.LQ"]], x))

# read in perditus locs
locs_per <- st_read(here("analysis", "output", "thinned_localities", "perditus_thinned.gpkg"))

# crop to perditus locs with buffer
## larger buffer than used to train SDM to allow for dispersal
mcp_per <- st_convex_hull(st_union(locs_per)) %>%
  st_buffer(dist = units::set_units(2.0, degree)) %>%
  terra::vect()

projections_per_crop <- crop(rast(projections_per), mcp_per) %>%
  mask(mcp_per)

# plot stability of environment through time
plot(stdev(projections_per_crop, na.rm = TRUE))
points(vect(locs_per))

# write predictions to file
writeRaster(projections_per_crop, here("analysis", "output", "sdm_projections", "projections_per.tif"))





