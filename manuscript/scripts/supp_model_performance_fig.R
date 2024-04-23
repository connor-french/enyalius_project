library(tidyverse)
library(pROC)
library(glue)
library(ggtext)
library(here)
library(patchwork)


# 0. set plotting theme ---------------------------------------------------

theme_set(theme_bw() + theme(
  plot.title = element_markdown(size = 16),
  axis.title = element_text(size = 12)
))

# PARAMETER ESTIMATION ---------------------------------------------------------------


# 1. load data ------------------------------------------------------------

# function to read in model results for parameter estimation
## test predictions begin with regression_test_preds_ and observed values begin with y_test_
read_model_results_reg <- function(species) {

  # if the species is cat or per, the number of admixture pops is 3, otherwise it is 2
  if (species == "cat" | species == "per") {
    num_pops <- 3
  } else {
    num_pops <- 2
  }

  admix_n_names <- paste0("admix_n_a", 1:num_pops)

  test_obs <-
    read_csv(
      here(
        "analysis",
        "output",
        "parameter_estimation",
        glue("reg_tuning_{species}"),
        glue("y_test_{species}.csv")
      )
    )

  test_pred <-
    read_csv(
      here(
        "analysis",
        "output",
        "parameter_estimation",
        glue("reg_tuning_{species}"),
        glue("test_preds.csv")
      ),
      skip = 0,
      col_names = c("max_k", "Nm", admix_n_names)
    )

  return(list(test_obs = test_obs, test_pred = test_pred))
}

pred_reg_cat <- read_model_results_reg("cat")
pred_reg_ihe <- read_model_results_reg("ihe")


# 5. observed vs predicted plots --------------------------------------------------------

plot_obs_vs_pred <- function(test_obs, test_pred, species, var) {
  if (species == "ihe") {
    sp_full <-  "E. iheringii"
    fill_col <- "#7CCE00"
  } else if (species == "cat") {
    sp_full <- "E. catenatus"
    fill_col <- "#CCA200"
  } else {
    sp_full <- "E. perditus"
    fill_col <- "#013220"
  }

  obs_pred <- tibble(
    obs = test_obs,
    pred = test_pred
  )

  obs_vs_pred <- ggplot(data = obs_pred, aes(x = pred, y = obs)) +
    geom_point(size = 2, fill = fill_col, shape = 21, alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      title = glue("{var}"),
      x = "Estimated",
      y = "Observed"
    )
  return(obs_vs_pred)
}

op_cat_a1 <-
  plot_obs_vs_pred(pred_reg_cat$test_obs$admix_n_a1,
                   pred_reg_cat$test_pred$admix_n_a1,
                   "cat",
                   "Admix N a1")
op_cat_a2 <-
  plot_obs_vs_pred(pred_reg_cat$test_obs$admix_n_a2,
                   pred_reg_cat$test_pred$admix_n_a2,
                   "cat",
                   "Admix N a2")

op_cat_a3 <-
  plot_obs_vs_pred(pred_reg_cat$test_obs$admix_n_a3,
                   pred_reg_cat$test_pred$admix_n_a3,
                   "cat",
                   "Admix N a3")



op_ihe_a1 <-
  plot_obs_vs_pred(
    test_obs = pred_reg_ihe$test_obs$admix_n_a1,
    test_pred = pred_reg_ihe$test_pred$admix_n_a1,
    "ihe",
    var = "Admix N a1"
  )

op_ihe_a2 <-
  plot_obs_vs_pred(
    test_obs = pred_reg_ihe$test_obs$admix_n_a2,
    test_pred = pred_reg_ihe$test_pred$admix_n_a2,
    "ihe",
    var = "Admix N a2"
  )

# COMBINE PLOTS  ---------------------------------------------------

layout <- "
AABBCC
#EEFF#
"

combo_plot <-
  op_cat_a1 + op_cat_a2 + op_cat_a3 + op_ihe_a1 + op_ihe_a2 +
  plot_layout(design = layout) +
  plot_annotation(
    title = "**Model performance**",
    tag_levels = "A",
    tag_suffix = ")",
    theme = theme(plot.title = element_markdown(size = 20))
  )

ggsave(here("manuscript", "figures", "supp_model_performance_fig.png"),
       combo_plot,
       width = 30,
       height = 15,
       units = "cm")

