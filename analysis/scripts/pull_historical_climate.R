# download and crop historical rasters from CHELSA
library(terra)
library(dplyr)
pull_historical_climate <- function(in_file, out_file, resample_template, crop_mcp) {
  r <- terra::rast(in_file)

  v <- terra::resample(r, resample_template)

  c <- terra::crop(v, crop_mcp) %>%
    terra::mask(crop_mcp)

  writeRaster(c, out_file, overwrite=TRUE)

}


