# script to generate the main SDM figure for the manuscript
# perditus is removed for dissertation

# 1. load packages ------------------------------------------------------------
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(here)
library(patchwork)
library(rnaturalearth)


# 2. read in localities ------------------------------------------------------
## sdm_locs prefix indicates localities used for SDM, after thinning
## snmf prefix indicates population assignment from sNMF

gen_locs <- read_csv(here("analysis", "data", "enyalius_locs_genetics.csv"))

### E. iheringii

#### read in snmf locs
snmf_ihe <- read_csv(here("analysis", "output", "snmf", "iheringii", "pop_assignment_a10_k2.csv")) %>%
  filter(!duplicated(id)) %>%
  select(id, pop = likely_assignment)

#### read in sdm locs
sdm_locs_ihe <- read_sf(here("analysis", "output", "thinned_localities", "iheringii_thinned.gpkg")) %>%
  left_join(snmf_ihe, by = c("id_code" = "id")) %>%
  # one sdm individual was not assigned to a population, but there is a representative individual in the snmf data
  mutate(
    pop = if_else(
      str_detect(id_code, "saop"), "pop_2", pop
    ),
    is_seq = if_else(!is.na(pop), "Yes", "No")
  )

### E. catenatus
gen_locs_cat <- gen_locs %>%
  filter(species == "catenatus") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

sdm_locs_cat <- read_sf(here("analysis", "output", "thinned_localities", "catenatus_thinned.gpkg"))

snmf_cat <- read_csv(here("analysis", "output", "snmf", "catenatus", "pop_assignment_a10.csv")) %>%
  filter(!duplicated(id)) %>%
  select(id, pop = likely_assignment)


sdm_locs_cat <- read_sf(here("analysis", "output", "thinned_localities", "catenatus_thinned.gpkg")) %>%
  left_join(snmf_cat, by = c("id_code" = "id")) %>%
  mutate(
    is_seq = if_else(!is.na(pop), "Yes", "No")
  )

### E. perditus
# gen_locs_per <- gen_locs %>%
#   filter(species == "perditus") %>%
#   st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
#
# sdm_locs_per <- read_sf(here("analysis", "output", "thinned_localities", "perditus_thinned.gpkg"))



# 3. read in SDM projections ----------------------------------------------

proj_ihe <- rast(here("analysis", "output", "sdm_projections", "projections_ihe.tif"))

proj_cat <- rast(here("analysis", "output", "sdm_projections", "projections_cat.tif"))

## in case I also display the catenatus projections with the inland removed
proj_cat_rm <- rast(here("analysis", "output", "sdm_projections", "projections_cat_inland-removed.tif"))

# proj_per <- rast(here("analysis", "output", "sdm_projections", "projections_per.tif"))

# 4. make individual plots --------------------------------------------------------
# present-day climate
## raw suitabilities
plot_sdm <- function(species, sdm, locs, sp_color, plot_type = "raw", leg_pos = "left") {

  if (plot_type == "raw") {
    sdm_plot <- ggplot() +
      geom_spatraster(data = sdm) +
      geom_sf(data = locs, fill = sp_color, aes(shape = pop), size = 2) +
      scale_fill_whitebox_c("muted") +
      scale_shape_manual(values = c(24, 21, 22), na.value = 13, guide = "none") +
      labs(title = species, fill = "Suitability") +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
      theme_bw() +
      theme(
        legend.position = leg_pos,
        legend.key.width = unit(1.0, 'cm'),
        legend.title = element_text(size = 12),
        plot.title = element_text(face = "bold.italic", size = 18)
      )
  } else if (plot_type == "thresh") {
    min_sdm <- terra::extract(sdm, locs, raw = TRUE) %>%
      min()

    sdm[sdm >= min_sdm] <- 1L
    sdm[sdm < min_sdm] <- 0

    sdm_plot <- ggplot() +
      geom_spatraster(data = sdm) +
      geom_sf(data = locs, fill = sp_color, aes(shape = pop), size = 2) +
      scale_fill_whitebox_c("muted") +
      scale_shape_manual(values = c(24, 21, 22), na.value = 13, guide = "none") +
      labs(title = species) +
      theme_bw() +
      theme(
        legend.position = leg_pos,
        legend.key.width = unit(0.5, 'cm'),
        legend.title = element_blank(),
        title = element_text(face = "italic")
      )
  } else {
    stop("plot_type must be either 'raw' or 'thresh'")
  }

  return(sdm_plot)

}

raw_ihe <-
  plot_sdm(
    species = "E. iheringii",
    sdm = proj_ihe$y0k,
    locs = sdm_locs_ihe,
    sp_color = "#0A5DA2",
    leg_pos = "bottom",
    plot_type = "raw"
  )

raw_cat <-
  plot_sdm(
    species = "E. catenatus",
    sdm = proj_cat_rm$`projections_cat_inland-removed_1`,
    locs = sdm_locs_cat %>% filter(longitude > -40.3),
    sp_color = "#DE8E07",
    leg_pos = "none",
    plot_type = "raw"
  )

#raw_per <- plot_sdm(species = "E. perditus", sdm = proj_per$y0k, locs = sdm_locs_per)


thresh_ihe <-
  plot_sdm(
    species = "",
    sdm = proj_ihe$y0k,
    locs = sdm_locs_ihe,
    sp_color = "#0A5DA2",
    plot_type = "thresh",
    leg_pos = "none"
  )

thresh_cat <- plot_sdm(
  species = "",
  proj_cat_rm$`projections_cat_inland-removed_1`,
  locs = sdm_locs_cat %>% filter(longitude > -40.3),
  sp_color = "#DE8E07",
  plot_type = "thresh",
  leg_pos = "none"
)

#thresh_per <- plot_sdm(species = "", sdm = proj_per$y0k, locs = sdm_locs_per)


# 5. make a geographic map for species range orientation based on their sdm projections  ------------------

# function to create an sf polygon outline of the sdm projections



plot_range <- function() {
  brazil <- ne_countries(scale = "medium", returnclass = "sf", country = "Brazil")

  range_ihe <- terra::mask(proj_ihe$y0k, proj_ihe$y0k > 0) %>%
    terra::as.polygons() %>%
    terra::aggregate()

  range_cat <- terra::mask(proj_cat_rm$`projections_cat_inland-removed_1`, proj_cat_rm$`projections_cat_inland-removed_1` > 0) %>%
    terra::as.polygons() %>%
    terra::aggregate()

  # range_per <- terra::mask(proj_per$y0k, proj_per$y0k > 0) %>%
  #   terra::as.polygons() %>%
  #   terra::aggregate()

  range_plot <- ggplot() +
    geom_sf(data = brazil) +
    geom_sf_text(
      data = brazil,
      aes(label = "bold('Brazil')"),
      parse = TRUE,
      nudge_y = -7,
      nudge_x = -18,
      size = 8
    ) +
    geom_sf(data = range_ihe, fill = "#0A5DA2", alpha = 1) +
    geom_sf_text(
      data = range_ihe,
      aes(label = "italic('E. iheringii')"),
      parse = TRUE,
      nudge_y = -0.5,
      nudge_x = 10,
      size = 5
    ) +
    geom_sf(data = range_cat, fill = "#DE8E07", alpha = 1) +
    geom_sf_text(
      data = range_cat,
      aes(label = "italic('E. catenatus')"),
      parse = TRUE,
      nudge_y = -0.5,
      nudge_x = -11,
      size = 5
    ) +
    # geom_sf(data = range_per, fill = "green", alpha = 0.5) +
    # geom_sf_text(data = range_per, aes(label = "italic('E. perditus')"), parse = TRUE, nudge_y = -0.5, nudge_x = 8, size = 5) +
    theme_void() +
    theme(legend.position = "none")

  return(range_plot)
}

range_map <- plot_range()

# 6. combine plots --------------------------------------------------------
full_plot <- (raw_cat + thresh_cat) /
(raw_ihe + thresh_ihe) /
# (raw_per + thresh_per) +
  range_map +
  plot_layout(widths = c(2, 2, 1)) +
  plot_annotation(tag_levels = "A", tag_suffix = ")")


# i further modify some positioning in adobe illustrator
ggsave(
  "sdm-main-fig.svg",
  full_plot,
  width = 20,
  height = 29,
  units = "cm",
  bg = "transparent",
  path = here("manuscript", "figures")
)










