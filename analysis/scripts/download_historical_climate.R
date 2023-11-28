library(terra)
library(dplyr)
library(stringr)
library(sf)
library(here)
library(readr)

source(here("analysis", "scripts", "pull_historical_climate.R"))

# script I used to download historical climate layers from CHELSA
# uses the pull_historical_climate functions from the pull_historical_climate.R script

template <- rast(here("analysis", "data", "biomes", "S1biomes_0k.asc")) %>%
  aggregate(2)

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

bioclim_nums <-
  c("bio03")



for (i in bioclim_nums) {
  for (j in seq(-60, -30, 10)) {
    bioclim_file <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_", i, "_", j, "_V1.0.tif")
    bioclim_out <- paste0("/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/enyalius-phylogeography/enyalius/analysis/output/historical_climate/", i, "_", j, ".tif")

    pull_historical_climate(in_file = bioclim_file, out_file = bioclim_out, resample_template = template, crop_mcp = mcp_all)

  }
}


