# spatial coalescent summary table
library(tidyverse)
library(glue)
library(here)


# 1. read in demographic summaries ----------------------------------------

read_demo <- function(species) {
  demo_sum_lin <- list.files(here("analysis", "output", "simulations", glue("sims_{species}_linear")), pattern = "demo", full.names = TRUE) %>%
    map_df(read_csv) %>%
    mutate(model = "linear")

  demo_sum_thresh <- list.files(here("analysis", "output", "simulations", glue("sims_{species}_threshold")), pattern = "demo", full.names = TRUE) %>%
    map_df(read_csv) %>%
    mutate(model = "threshold")

  demo_sum <- bind_rows(
    demo_sum_lin,
    demo_sum_thresh
  )

  return(demo_sum)
}

## catenatus
demo_sum_cat <- read_demo("cat")

## iheringii
demo_sum_ihe <- read_demo("ihe")

## perditus


# 2. create summary table -------------------------------------------------
## catenatus
demo_sum_cat %>%
  group_by(model) %>%
  summarize(
    total_inds_curr_max = max(total_inds_curr),
    total_inds_curr_min = min(total_inds_curr),
    total_inds_largest_max = max(total_inds_largest),
    total_inds_largest_min = min(total_inds_largest),
    avg_deme_size_max = max(avg_deme_size_avg),
    avg_deme_size_min = min(avg_deme_size_avg),
    num_occ_average_max = max(num_occ_average),
    num_occ_average_min = min(num_occ_average)
  ) %>% View()

## iheringii
demo_sum_ihe %>%
  group_by(model) %>%
  summarize(
    total_inds_curr_max = max(total_inds_curr),
    total_inds_curr_min = min(total_inds_curr),
    total_inds_largest_max = max(total_inds_largest),
    total_inds_largest_min = min(total_inds_largest),
    avg_deme_size_max = max(avg_deme_size_avg),
    avg_deme_size_min = min(avg_deme_size_avg),
    num_occ_average_max = max(num_occ_average),
    num_occ_average_min = min(num_occ_average)
  ) %>% View()

