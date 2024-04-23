# maps of genetic diversity across the landscape
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(patchwork)
library(here)


# 1. load data ---------------------------------------------------------

## load rasters
gdrast_ihe <- rast(here("analysis", "output", "gendiv_rasters", "gendiv_ihe.tif"))
bg_ihe <- rast(here("analysis", "output", "sdm_projections", "projections_ihe.tif"))

gdrast_cat <- rast(here("analysis", "output", "gendiv_rasters", "gendiv_cat.tif"))
bg_cat <- rast(here("analysis", "output", "sdm_projections", "projections_cat_inland-removed.tif"))

## get outline of the study area
outline_ihe <- bg_ihe$y0k %>%
  as.polygons() %>%
  aggregate()

outline_cat <- bg_cat$`projections_cat_inland-removed_1` %>%
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


# 2. plot gendiv rasters ---------------------------------------------------------
plot_gendiv <- function(rast, outline, locs, species, sp_color) {
  if (species == "E. catenatus") {
    breaks <- c(0.0015, 0.0045)
  } else if (species == "E. iheringii") {
    breaks <- c(0.001, 0.0086)
  } else if (species == "E. perditus") {
    breaks <- c(0.001, 0.005)
  }

  ggplot() +
    geom_spatvector(data = outline, fill = "gray") +
    geom_spatraster(data = rast) +
    scale_fill_whitebox_c("viridi", breaks = breaks, labels = c("Low", "High")) +
    geom_sf(data = locs, fill = sp_color, size = 1, shape = 21) +
    guides(fill = guide_colorbar(title = "Genetic\ndiversity")) +
    labs(title = species) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold.italic"),
      axis.text.x = element_text(angle = 45, vjust = 0.75)
      )
}

gendiv_plot_ihe <-
  plot_gendiv(gdrast_ihe, outline_ihe, gen_locs_ihe, species = "E. iheringii", sp_color = "#7CCE00")

gendiv_plot_cat <- plot_gendiv(gdrast_cat, outline_cat, gen_locs_cat, species = "E. catenatus", sp_color = "#CCA200") + theme(legend.position = "none")

# 4. calculate stability ---------------------------------------------------------
calc_stab <- function(r, thresh) {

  ## apply the threshold for all layers in the stack
  r_thresh <- r > thresh

  ## sum the number of layers that are above the threshold
  r_sum <- sum(r_thresh)

  ## threshold so zero values are NA
  r_sum[r_sum == 0] <- NA

  return(r_sum)
}

thresh_ihe <- 0.23516189
thresh_cat <- 0.36185417

stab_ihe <- calc_stab(bg_ihe, thresh_ihe)
stab_cat <- calc_stab(bg_cat, thresh_cat)

# 4. plot stability ---------------------------------------------------------
plot_stab <- function(stab, outline, locs, species, sp_color) {
  ggplot() +
    geom_spatvector(data = outline, fill = "gray") +
    geom_spatraster(data = stab) +
    scale_fill_whitebox_c("viridi", breaks = c(4, 20), labels = c("Low", "High")) +
    geom_sf(data = locs, fill = sp_color, size = 1, shape = 21) +
    guides(fill = guide_colorbar(title = "Stability")) +
    labs(title = species) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 0.75)
    )
}

stab_plot_cat <- plot_stab(stab_cat, outline_cat, gen_locs_cat, species = "E. catenatus", sp_color = "#CCA200") + theme(legend.position = "none")

stab_plot_ihe <- plot_stab(stab_ihe, outline_ihe, gen_locs_ihe, species = "E. iheringii", sp_color = "#7CCE00")



# 5. full plot ------------------------------------------------------------

(gendiv_plot_cat +
  gendiv_plot_ihe) /
  (stab_plot_cat +
    stab_plot_ihe) +
  plot_annotation(tag_levels = "A", tag_suffix = ")")


# 4. save plot ------------------------------------------------------------
ggsave(here("manuscript", "figures", "gendiv_maps.png"), width = 15, height = 15, units = "cm", dpi = 300)




