# script to plot a rough assessment of climate stability
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(here)
library(patchwork)


# 1. load the data --------------------------------------------------------

proj_ihe <- rast(here("analysis", "output", "sdm_projections", "projections_ihe.tif"))
proj_cat <- rast(here("analysis", "output", "sdm_projections", "projections_cat_inland-removed.tif"))

## get outline of the study area
outline_ihe <- proj_ihe$y0k %>%
  as.polygons() %>%
  aggregate()

outline_cat <- proj_cat$`projections_cat_inland-removed_1` %>%
  as.polygons() %>%
  aggregate()

## read in localities
gen_locs <- read_csv(here("analysis", "data", "enyalius_locs_genetics.csv"))

### E. iheringii
gen_locs_ihe <- gen_locs %>%
  filter(species == "iheringii") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

### E. catenatus
inland_inds = c("cat_jeq_MTR17108", "cat_jeq_MTR17172", "cat_jeq_MTR17377", "cat_mig_MTR19889", "cat_mig_RPD188", "cat_rem_MTR19915", "cat_rem_RPD128")

gen_locs_cat <- gen_locs %>%
  filter(
    species == "catenatus",
    !id_code %in% inland_inds
  ) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)



# 2. apply a threshold and sum the layers --------------------------------------------------------
# iheringii
## threshold applied to the simulations
thresh_ihe <- 0.23516189

## apply the threshold for all layers in the stack
proj_ihe_thresh <- proj_ihe > thresh_ihe

## sum the number of layers that are above the threshold
proj_ihe_sum <- sum(proj_ihe_thresh)

## threshold so zero values are NA
proj_ihe_sum[proj_ihe_sum == 0] <- NA

# catenatus

## threshold applied to the simulations
thresh_cat <- 0.36185417

## apply the threshold for all layers in the stack
proj_cat_thresh <- proj_cat > thresh_cat

## sum the number of layers that are above the threshold
proj_cat_sum <- sum(proj_cat_thresh)

## threshold so zero values are NA
proj_cat_sum[proj_cat_sum == 0] <- NA

# 3. plot the maps --------------------------------------------------------

plot_stab <- function(rast, outline, locs, species, sp_color) {
  ggplot() +
    geom_spatvector(data = outline, fill = "gray") +
    geom_spatraster(data = rast) +
    scale_fill_whitebox_c("viridi") +
    geom_sf(data = locs, fill = sp_color, size = 1, shape = 21) +
    guides(fill = guide_colorbar(title = "Stability")) +
    labs(title = species) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold.italic"))
}

plot_stab(proj_cat_sum, outline_cat, gen_locs_cat, "E. catenatus", "#CCA200") +
  plot_stab(proj_ihe_sum, outline_ihe, gen_locs_ihe, "E. iheringii", "#7CCE00") +
  plot_annotation(tag_levels = "A", tag_suffix = ")") +
  plot_layout(guides = "collect")


# 4. save the plot --------------------------------------------------------
ggsave(here("manuscript", "figures", "stability_maps.png"), width = 10, height = 5, dpi = 300)


