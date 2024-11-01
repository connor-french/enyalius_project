# figure for the species range map and elevational distribution
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(rnaturalearth)
library(patchwork)
library(ggnewscale)
library(ggridges)
library(here)



# 1. load the data --------------------------------------------------------
## read in localities
## inland individuals to exclude
inland_inds = c("cat_jeq_MTR17108", "cat_jeq_MTR17172", "cat_jeq_MTR17377", "cat_mig_MTR19889", "cat_mig_RPD188", "cat_rem_MTR19915", "cat_rem_RPD128")

## read in the localities and filter for iheringii and catenatus
gen_locs <- read_csv(here("analysis", "data", "enyalius_locs_genetics.csv")) %>%
  filter(species == "iheringii" | species == "catenatus",
         !id_code %in% inland_inds) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


## read in DEM for extracting elevations
dem <- rast(here("manuscript", "data", "af_dem.tif"))

## read in coarser DEM for plotting (the finer resolution DEM leads to spurious white lines in the plots for some reason)
dem_bio <- rast(here("manuscript", "data", "wc2.1_2.5m_elev.tif")) %>%
  crop(ext(dem))

## read in the af boundary
af <- terra::vect(here("analysis", "data", "atlantic_forest", "atlantic_forest.geojson")) %>%
  aggregate()

# 2. extract elevation from points
gen_locs <- gen_locs %>%
  mutate(elevation = extract(dem, .)[,2])

# 3. plot the data --------------------------------------------------------
## map


sp_map <- ggplot() +
  geom_spatraster(data = dem_bio) +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent") +
  labs(fill = "Elev. (m)") +
  geom_spatvector(data = af, fill = "darkgreen", color = "darkgreen", linewidth = 0.7, alpha = 0.2) +
  geom_spatvector(data = af, fill = "transparent", color = "darkgreen", linewidth = 0.7) +
  new_scale_fill() +
  geom_spatvector(data = gen_locs, aes(fill = species), shape = 21, size = 2.5) +
  scale_fill_manual(values = c("#DE8E07", "#0A5DA2"), guide = "none") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.365, 0.15)
  )



## inset map
brazil <- ne_countries(scale = "medium", country = "brazil", returnclass = "sf")

brazil_inset <- ggplot() +
  geom_sf(data = brazil, fill = "white") +
  geom_sf_text(data = brazil, aes(label = name), size = 7, nudge_x = -2.85, nudge_y = 4.25) +
  geom_spatvector(data = af, fill = "darkgreen", color = "darkgreen") +
  theme_void()

## density plot
### specifying the beginning and end points of the curve to connect the label to the plot
splab_df <- tibble(
  species = c("E. iheringii", "E. catenatus"),
  x = c(10, 1000),
  y = c(1.65, 2.75)
)

sp_density <- gen_locs %>%
  mutate(species_num = if_else(species == "catenatus", 2, 1)) %>%
  ggplot() +
  geom_density_ridges(aes(x = elevation, y = species_num, group = species_num, fill = species), scale = 1) +
  scale_fill_manual(values = c("#DE8E07", "#0A5DA2"), guide="none") +
  geom_text(data = splab_df, aes(x = x, y = y, label = species), fontface = "bold.italic", size = 4.3) +
  geom_curve(x = 660, y = 2.75, xend = 380, yend = 2.55, curvature = 0.2, color = "black") +
  geom_curve(x = 340, y = 1.65, xend = 600, yend = 1.55, curvature = -0.2, color = "black") +
  labs(x = "Elevation (m)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12)
  )

sp_combo <- sp_map +
  inset_element(sp_density, left = 0.43, bottom = 0, right = 1, top = 0.38) +
  inset_element(brazil_inset, left = 0.05, bottom = 0.55, right = 0.45, top = 1.05, clip = FALSE)


# 4. save the figure ------------------------------------------------------
ggsave(here("manuscript", "figures", "range_map.png"), sp_combo, width = 30, height = 20, units = "cm", dpi = 600)

