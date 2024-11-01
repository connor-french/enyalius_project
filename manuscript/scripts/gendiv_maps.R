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
snmf_ihe <- read_csv(here("analysis", "output", "snmf", "iheringii", "pop_assignment_a10_k2.csv")) %>%
  filter(!duplicated(id)) %>%
  select(id, pop = likely_assignment)

gen_locs_ihe <- gen_locs %>%
  filter(species == "iheringii") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  left_join(snmf_ihe, by = c("id_code" = "id"))

### E. catenatus
inland_inds = c("cat_jeq_MTR17108", "cat_jeq_MTR17172", "cat_jeq_MTR17377", "cat_mig_MTR19889", "cat_mig_RPD188", "cat_rem_MTR19915", "cat_rem_RPD128")

#### get snmf pop assignment
snmf_cat <- read_csv(here("analysis", "output", "snmf", "catenatus", "pop_assignment_a10.csv")) %>%
  filter(!duplicated(id)) %>%
  select(id, pop = likely_assignment)

gen_locs_cat <- gen_locs %>%
  filter(
    species == "catenatus",
    !id_code %in% inland_inds
         ) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  left_join(snmf_cat, by = c("id_code" = "id"))



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
    geom_sf(data = locs, fill = sp_color, aes(shape = pop), size = 2) +
    scale_shape_manual(values = c(24, 21, 22), na.value = 13, guide = "none") +
    guides(fill = guide_colorbar(title = "Genetic\ndiversity")) +
    labs(title = species) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold.italic"),
      axis.text.x = element_text(angle = 45, vjust = 0.75)
      )
}

gendiv_plot_ihe <-
  plot_gendiv(gdrast_ihe, outline_ihe, gen_locs_ihe, species = "E. iheringii", sp_color = "#0A5DA2")

gendiv_plot_cat <- plot_gendiv(gdrast_cat, outline_cat, gen_locs_cat, species = "E. catenatus", sp_color = "#DE8E07") + theme(legend.position = "none")

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
    geom_sf(data = locs, fill = sp_color, aes(shape = pop), size = 2) +
    scale_shape_manual(values = c(24, 21, 22), na.value = 13, guide = "none") +
    guides(fill = guide_colorbar(title = "Stability")) +
    labs(title = species) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 0.75)
    )
}

stab_plot_cat <- plot_stab(stab_cat, outline_cat, gen_locs_cat, species = "E. catenatus", sp_color = "#DE8E07") + theme(legend.position = "none")

stab_plot_ihe <- plot_stab(stab_ihe, outline_ihe, gen_locs_ihe, species = "E. iheringii", sp_color = "#0A5DA2")


# 5. landscape size over time ------------------------------------------------------------
# projections
# posterior N
set.seed(4222)
n_post_cat <- runif(1000, 1162, 4970)
n_post_ihe <- runif(1000, 275, 2987)

# point estimate
n_point_cat <- 3083
n_point_ihe <- 1588

# get the sum of the landscape-wide sdm values after multiplying them by the posterior N for each timestep

thresh_ihe <- 0.23516189
thresh_cat <- 0.36185417

# sum the values of the demes for each layer
sum_layer <- function(l) {

  return(sum(values(l), na.rm = TRUE))
}

sum_vals <- function(p, post, thresh) {

  p[p < thresh] <- 0L
  p[p >= thresh] <- 1L

  p <- p * post
  s_ts <- map_dbl(1:nlyr(p), \(x) sum_layer(p[[x]])) %>%
    as_tibble_col(column_name = paste0("n_", round(post, 0)))

  return(s_ts)
}

sum_demes_cat <- map(n_post_cat, sum_vals, p = bg_cat, thresh = thresh_cat) %>%
  list_cbind() %>%
  mutate(time = 0:22) %>%
  pivot_longer(cols = starts_with("n_"), names_to = "post_n", values_to = "sum_n")

sum_demes_ihe <- map(n_post_ihe, sum_vals, p = bg_ihe, thresh = thresh_ihe) %>%
  list_cbind() %>%
  mutate(time = 0:22) %>%
  pivot_longer(cols = starts_with("n_"), names_to = "post_n", values_to = "sum_n")


# create dataframes for the point estimates
point_cat <- tibble(time = 0:22, post_n = "point", sum_n = pluck(sum_vals(bg_cat, n_point_cat, thresh = thresh_cat), 1))
point_ihe <- tibble(time = 0:22, post_n = "point", sum_n = pluck(sum_vals(bg_ihe, n_point_ihe, thresh = thresh_ihe), 1))


# plot cat, with the point estimate highlighted
ts_plot_cat <- ggplot() +
  geom_line(data = sum_demes_cat, aes(x = time, y = sum_n, group = post_n), color = "gray30") +
  geom_line(data = point_cat, aes(x = time, y = sum_n), color = "red", linewidth = 2) +
  labs(x = "Time (kya)", y = "Total N") +
  scale_y_continuous(labels = scales::label_scientific(), limits = c(0, 9e06)) +
  theme_bw()

# plot ihe, with the point estimate highlighted
ts_plot_ihe <- ggplot() +
  geom_line(data = sum_demes_ihe, aes(x = time, y = sum_n, group = post_n), color = "gray30") +
  geom_line(data = point_ihe, aes(x = time, y = sum_n), color = "red", linewidth = 2) +
  labs(x = "Time (kya)", y = "Total N") +
  scale_y_continuous(labels = scales::label_scientific(), limits = c(0, 9e06)) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())


# 6. full plot ------------------------------------------------------------

(gendiv_plot_cat +
  gendiv_plot_ihe) /
  (stab_plot_cat +
    stab_plot_ihe) /
  (ts_plot_cat +
    ts_plot_ihe) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") +
  plot_layout(heights = c(2, 2, 1))


# 4. save plot ------------------------------------------------------------
ggsave(here("manuscript", "figures", "gendiv_maps.svg"), width = 15, height = 20, units = "cm", dpi = 300)




