
# load packages -----------------------------------------------------------

# PCA exploration of simulation results
library(tidyverse)
library(tidymodels)
library(applicable)
library(glue)
library(here)
library(patchwork)


# 1. define functions -----------------------------------------------------

# three letter code for species
read_sims <- function(species) {
  params_lin <-
    list.files(
      here("analysis", "output", "simulations", glue("sims_{species}_linear")),
      pattern = "params.csv",
      full.names = TRUE
    ) %>%
    map_df(read_csv, show_col_types = FALSE)

  params_thresh <-
    list.files(
      here("analysis", "output", "simulations", glue("sims_{species}_threshold")),
      pattern = "params.csv",
      full.names = TRUE
    ) %>%
    map_df(read_csv, show_col_types = FALSE)

  demo_sum_lin <-
    list.files(
      here("analysis", "output", "simulations", glue("sims_{species}_linear")),
      pattern = "demo_summary.csv",
      full.names = TRUE
    ) %>%
    map_df(read_csv, show_col_types = FALSE)

  demo_sum_thresh <-
    list.files(
      here("analysis", "output", "simulations", glue("sims_{species}_threshold")),
      pattern = "demo_summary.csv",
      full.names = TRUE
    ) %>%
    map_df(read_csv, show_col_types = FALSE)

  sumstats_lin <-
    list.files(
      here("analysis", "output", "simulations", glue("sims_{species}_linear")),
      pattern = "sumstats.csv",
      full.names = TRUE
    ) %>%
    map_df(read_csv, show_col_types = FALSE) %>%
    select(
      -morans_i,
      -starts_with("pi_pop"),
      -contains("h3"),
      -param_id,
      -sim_id
    )

  sumstats_thresh <-
    list.files(
      here("analysis", "output", "simulations", glue("sims_{species}_threshold")),
      pattern = "sumstats.csv",
      full.names = TRUE
    ) %>%
    map_df(read_csv, show_col_types = FALSE) %>%
    select(
      -morans_i,
      -starts_with("pi_pop"),
      -contains("h3"),
      -param_id,
      -sim_id
    )

  # full dataframes
  full_lin <- bind_cols(params_lin, sumstats_lin, demo_sum_lin) %>%
    mutate(model = "linear")

  full_thresh <- bind_cols(params_thresh, sumstats_thresh, demo_sum_thresh) %>%
    mutate(model = "threshold")

  full <- bind_rows(full_lin, full_thresh)

  # only sumstats
  full_ss <- bind_rows(sumstats_lin %>% mutate(model = "linear"), sumstats_thresh %>% mutate(model = "threshold"))



  return(list(full = full, sumstats = full_ss))

}

# read in and filter empirical summary statistics
read_empirical <- function(path){
  emp <- read_csv(path, show_col_types = FALSE) %>%
    select(-contains("fst"),
           -...1,
           -contains("h3"),
           -num_var) %>%
    rename_with(~str_remove(.x, "_dxy"), starts_with("ibd"))

  return(emp)
}



# 3. read in data ---------------------------------------------------------

## catenatus
ss_cat <- read_sims("cat")

emp_cat <- read_empirical(here("analysis", "output", "empirical_sumstats", "catenatus_sumstats.csv"))

## iheringii
ss_ihe <- read_sims("ihe")

emp_ihe <- read_empirical(here("analysis", "output", "empirical_sumstats", "iheringii_sumstats.csv"))


# 4. assess model extrapolation -------------------------------------------

## catenatus
# set up a tidymodels recipe
recipe_cat <-
  recipe(formula = model ~ ., data = ss_cat$sumstats) %>%
  step_naomit(starts_with("taj_d")) %>%
  step_zv(all_predictors())

# perform apd PCA
pca_cat <- apd_pca(recipe_cat, ss_cat$sumstats)

# get coverage score of empirical data
pca_score_cat <- score(pca_cat, emp_cat)

pca_score_cat %>% select(matches("PC00[1-2]"), contains("distance")) %>% pull(distance_pctl)

## iheringii
# set up a tidymodels recipe
recipe_ihe <-
  recipe(formula = model ~ ., data = ss_ihe$sumstats) %>%
  step_zv(all_predictors())

# perform PCA
pca_ihe <- apd_pca(recipe_ihe, ss_ihe$sumstats)

# get coverage score of empirical data
pca_score_ihe <- score(pca_ihe, emp_ihe)

pca_score_ihe %>% select(matches("PC00[1-2]"), contains("distance")) %>% pull(distance_pctl)


# 5. perform PCA for plotting ---------------------------------------------

## catenatus
pca_cat <- prcomp(ss_cat$sumstats %>% na.omit() %>% select(-model), center = TRUE, scale. = TRUE)

pca_df_cat <- ss_cat$sumstats %>%
  na.omit() %>%
  select(model) %>%
  bind_cols(pca_cat$x)

emp_pca <- predict(pca_cat, emp_cat) %>%
  as_tibble()

pca_plot_cat <- ggplot() +
  geom_point(data = pca_df_cat, aes(x = PC1, y = PC2), color = "gray") +
  geom_point(data = emp_pca, aes(x = PC1, y = PC2), color = "red")  +
  labs(title = "E. catenatus") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "italic")
    )


## iheringii
pca_ihe <- prcomp(ss_ihe$sumstats %>% na.omit() %>% select(-model), center = TRUE, scale. = TRUE)

pca_df_ihe <- ss_ihe$sumstats %>%
  na.omit() %>%
  select(model) %>%
  bind_cols(pca_ihe$x)

emp_pca <- predict(pca_ihe, emp_ihe) %>%
  as_tibble()

pca_plot_ihe <- ggplot() +
  geom_point(data = pca_df_ihe, aes(x = PC1, y = PC2), color = "gray") +
  geom_point(data = emp_pca, aes(x = PC1, y = PC2), color = "red")  +
  labs(title = "E. iheringii") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "italic")
  )

pca_plot_total <-
  pca_plot_cat + pca_plot_ihe

ggsave(here("manuscript", "figures", "gof-pca-fig.png"),
       pca_plot_total,
       width = 20,
       height = 10,
       units = "cm")




