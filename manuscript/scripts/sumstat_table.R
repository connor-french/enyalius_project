# summary statistic table for manuscript
library(tidyverse)
library(here)


## catenatus
ss_cat <- read_csv(here("analysis", "output", "empirical_sumstats", "catenatus_sumstats.csv")) %>%
  select(
    -contains("fst"),
    -...1,
    -num_var,
    -contains("h3"),
    -starts_with("pi_pop")
  ) %>%
  rename_with(~str_remove(.x, "_dxy"), starts_with("ibd"))


dxy_cat <- ss_cat %>%
  pivot_longer(any_of(starts_with("dxy")), names_to = "dxy_name", values_to = "dxy") %>%
  summarize(mean_dxy = mean(dxy))


ss_table_cat <- ss_cat %>%
  select(-starts_with("dxy")) %>%
  mutate(mean_dxy = dxy_cat$mean_dxy)

# write to clipboard for pasting into gdocs
clipr::write_clip(ss_table_cat)

## iheringii
ss_ihe <- read_csv(here("analysis", "output", "empirical_sumstats", "iheringii_sumstats.csv")) %>%
  select(
    -contains("fst"),
    -...1,
    -num_var,
    -contains("h3"),
    -starts_with("pi_pop")
  ) %>%
  rename_with(~str_remove(.x, "_dxy"), starts_with("ibd"))


dxy_ihe <- ss_ihe %>%
  pivot_longer(any_of(starts_with("dxy")), names_to = "dxy_name", values_to = "dxy") %>%
  summarize(mean_dxy = mean(dxy))


ss_table_ihe <- ss_ihe %>%
  select(-starts_with("dxy")) %>%
  mutate(mean_dxy = dxy_ihe$mean_dxy)

## perditus





