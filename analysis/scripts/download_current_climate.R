# script I used to download current climate layers from CHELSA
library(terra)
library(sf)
library(here)
library(readr)
library(dplyr)
library(stringr)

# download and crop historical rasters from CHELSA
pull_current_climate <- function(in_file, out_file, resample_template, crop_mcp) {

  r <- terra::rast(in_file)

  v <- terra::resample(r, resample_template)

  c <- terra::crop(v, crop_mcp) %>%
    terra::mask(crop_mcp)

  template_crop <- terra::crop(resample_template, crop_mcp) %>%
    terra::mask(crop_mcp)

  c <- c * template_crop

  writeRaster(c, out_file)

}


template <- rast(here("analysis", "data", "biomes", "S1biomes_0k.asc")) %>%
  aggregate(2)

# create a mask for restricting the bioclim coast to be the same as the biome template file
template[template > 0] <- 1
template[template < 1] <- 0

locs_df <- read_csv(here("analysis", "data", "enyalius_locs.csv")) %>%
  # the variable names are messy
  janitor::clean_names() %>%
  filter(!is.na(latitude), !is.na(longitude),
         # remove localities that are only mapped to a Google Earth municipality. The reserve in the "GPS of the reserve" locality is barely a km, so that is acceptable error relative to the environmental resolution
         #source_latlong != "Google Earth municipality"
  ) %>%
  # there are some lat/longs with spaces. Need to fix this
  mutate(longitude = str_remove(longitude, "\\s"),
         latitude = str_remove(latitude, "\\s")) %>%
  # convert to numeric
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude))

# convert to sf for spatial operations
locs_sf <- st_as_sf(locs_df,
                    coords = c("longitude", "latitude"),
                    crs = 4326)


mcp_all <- st_convex_hull(st_union(locs_sf)) %>%
  st_buffer(dist = units::set_units(2, degree)) %>%
  terra::vect()

# if the internet gets interrupted, you may need to re-run this loop, starting with the raster that failed to downloaded
bioclim_nums <- paste0("bio", 1:19)

for (i in bioclim_nums) {
    bioclim_file <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_", i, "_1981-2010_V.2.1.tif")

    bioclim_out <- file.path(here("analysis", "data", "current_climate_chelsa"), paste0(i,".tif"))

    pull_current_climate(in_file = bioclim_file, out_file = bioclim_out, resample_template = template, crop_mcp = mcp_all)

}
