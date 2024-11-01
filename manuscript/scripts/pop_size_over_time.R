# script to estimate the population size change over time and plot it
library(tidyverse)
library(terra)
library(here)

# projections
projections_cat <- rast(here("analysis", "output", "sdm_projections", "projections_cat_inland-removed.tif"))
projections_ihe <- rast(here("analysis", "output", "sdm_projections", "projections_ihe.tif"))


# posterior N
set.seed(4222)
n_post_cat <- runif(1000, 1162, 4970)
n_post_ihe <- runif(1000, 275, 2987)

# point estimate
n_point_cat <- 3083
n_point_ihe <- 1588

# get the sum of the landscape-wide sdm values after multiplying them by the posterior N for each timestep

# sum the values of the demes for each layer
sum_layer <- function(l) {
  return(sum(values(l), na.rm = TRUE))
}

sum_vals <- function(p, post) {
  p <- p * post
  s_ts <- map_dbl(1:nlyr(p), \(x) sum_layer(p[[x]])) %>%
    as_tibble_col(column_name = paste0("n_", round(post, 0)))

  return(s_ts)
}

sum_demes_cat <- map(n_post_cat, sum_vals, p = projections_cat) %>%
  list_cbind() %>%
  mutate(time = 0:22) %>%
  pivot_longer(cols = starts_with("n_"), names_to = "post_n", values_to = "sum_n")

sum_demes_ihe <- map(n_post_ihe, sum_vals, p = projections_ihe) %>%
  list_cbind() %>%
  mutate(time = 0:22) %>%
  pivot_longer(cols = starts_with("n_"), names_to = "post_n", values_to = "sum_n")


# create dataframes for the point estimates
point_cat <- tibble(time = 0:22, post_n = "point", sum_n = pluck(sum_vals(projections_cat, n_point_cat), 1))
point_ihe <- tibble(time = 0:22, post_n = "point", sum_n = pluck(sum_vals(projections_ihe, n_point_ihe), 1))


# plot cat, with the point estimate highlighted
ggplot() +
  geom_line(data = sum_demes_cat, aes(x = time, y = sum_n, group = post_n), color = "gray30") +
  geom_line(data = point_cat, aes(x = time, y = sum_n), color = "red", linewidth = 2) +
  labs(x = "Time (kya)", y = "Total N across the landscape") +
  lims(y = c(0, 9e06)) +
  theme_minimal()

# plot ihe, with the point estimate highlighted
ggplot() +
  geom_line(data = sum_demes_ihe, aes(x = time, y = sum_n, group = post_n), color = "gray30") +
  geom_line(data = point_ihe, aes(x = time, y = sum_n), color = "red", linewidth = 2) +
  labs(x = "Time (kya)", y = "Total N across the landscape") +
  lims(y = c(0, 9e06)) +
  theme_bw()

