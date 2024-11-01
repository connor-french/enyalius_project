# plot model performance results
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

# CLASSIFICATION ----------------------------------------------------------

# 1. load data ------------------------------------------------------------

# function to read in model results for classification
read_model_results <- function(species) {
  test_obs <-
    read_csv(
      here(
        "analysis",
        "output",
        "model_classification",
        glue("results_{species}"),
        glue("y_test_{species}.csv")
      ),
      col_names = "transformation"
    ) %>%
    mutate(transformation = as.integer(transformation))

  test_pred_bin <-
    read_csv(
      here(
        "analysis",
        "output",
        "model_classification",
        glue("results_{species}"),
        glue("test_pred_bin_{species}.csv")
      ),
      col_names = "transformation"
    ) %>%
    mutate(transformation = as.integer(transformation))

  test_pred_prob <-
    read_csv(
      here(
        "analysis",
        "output",
        "model_classification",
        glue("results_{species}"),
        glue("test_pred_prob_{species}.csv")
      ),
      col_names = c("prob_0", "prob_1")
    )

  return(list(test_obs = test_obs, test_pred_bin = test_pred_bin, test_pred_prob = test_pred_prob))
}

# read in the data
data_cat <- read_model_results("cat")
data_ihe <- read_model_results("ihe")


# 2. ROC curve ------------------------------------------------------------


# turn all of the code above into a function
plot_roc_curve <- function(test_obs, test_pred_prob, species) {
  roc_obj <- roc(test_obs$transformation, test_pred_prob$prob_0)
  auc_score <- auc(roc_obj)

  roc_df <- tibble(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  if (species == "ihe") {
    lincol <- "#0A5DA2"
    sp_full <-  "E. iheringii"
  } else if (species == "cat") {
    lincol <-  "#DE8E07"
    sp_full <- "E. catenatus"
  } else {
    lincol <-  "black"
  }

  auc_label_df <- tibble(
    x = 0.75,
    y = 0.35,
    lab = paste0("AUC: ", round(auc_score, 2))
  )

  roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(linewidth = 1.2, color = lincol) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    # inset label with the AUC score
    geom_label(data = auc_label_df,
     aes(x = x, y = y, label = lab),
      fill = "white",
      color = "black",
      size = 4,
      label.size = 0.5
    ) +
    labs(
      title = glue("*{sp_full}*"),
      x = "False Positive Rate",
      y = "True Positive Rate"
    )


  return(roc_plot)
}

roc_cat <- plot_roc_curve(data_cat$test_obs, data_cat$test_pred_prob, "cat")
roc_ihe <- plot_roc_curve(data_ihe$test_obs, data_ihe$test_pred_prob, "ihe")

# 3. confusion matrix -----------------------------------------------------
plot_confusion_matrix <- function(test_obs, test_pred_bin, species) {
  conf_mat <- table(test_obs$transformation, test_pred_bin$transformation)
  # make a confusion matrix of percentages
  conf_mat <- prop.table(conf_mat, margin = 1) * 100

  if (species == "ihe") {
    sp_full <-  "E. iheringii"
    fill_col <- "#0A5DA2"
  } else if (species == "cat") {
    sp_full <- "E. catenatus"
    fill_col <- "#DE8E07"
  } else {
    sp_full <- "species"
    fill_col <- "darkgreen"
  }
  conf_mat_plot <- ggplot(data = as.data.frame(conf_mat), aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = Freq), color = "black") +
    geom_text(aes(label = paste0(round(Freq, 2), "%"))) +
    scale_fill_gradient(low = "white", high = fill_col) +
    labs(
      title = glue("*{sp_full}*"),
      x = "Observed",
      y = "Estimated"
    ) +
    # change axis tick labels that are 0 to Linear and 1 to Threshold
    scale_x_discrete(labels = c("Linear", "Threshold")) +
    scale_y_discrete(labels = c("Linear", "Threshold")) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      legend.position = "none"
      )
  return(conf_mat_plot)
}

conf_mat_cat <- plot_confusion_matrix(data_cat$test_obs, data_cat$test_pred_bin, "cat")
conf_mat_ihe <- plot_confusion_matrix(data_ihe$test_obs, data_ihe$test_pred_bin, "ihe")

# PARAMETER ESTIMATION ---------------------------------------------------------------


# 4. load data ------------------------------------------------------------

# function to read in model results for parameter estimation
## test predictions begin with regression_test_preds_ and observed values begin with y_test_
read_model_results_reg <- function(species) {

  # if the species is cat or per, the number of admixture pops is 3, otherwise it is 2
  if (species == "cat" | species == "per") {
    num_pops <- 3
  } else {
    num_pops <- 2
  }

  admix_n_names <- paste0("admix_n_", 1:num_pops)

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
    fill_col <- "#0A5DA2"
  } else if (species == "cat") {
    sp_full <- "E. catenatus"
    fill_col <- "#DE8E07"
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

op_cat_maxk <-
  plot_obs_vs_pred(pred_reg_cat$test_obs$max_k,
                   pred_reg_cat$test_pred$max_k,
                   "cat",
                   "Max deme size")
op_cat_nm <-
  plot_obs_vs_pred(pred_reg_cat$test_obs$Nm,
                   pred_reg_cat$test_pred$Nm,
                   "cat",
                   "*Nm*")

op_ihe_maxk <-
  plot_obs_vs_pred(
    test_obs = pred_reg_ihe$test_obs$max_k,
    test_pred = pred_reg_ihe$test_pred$max_k,
    "ihe",
    var = "Max deme size"
  )

op_ihe_nm <-
  plot_obs_vs_pred(
    test_obs = pred_reg_ihe$test_obs$Nm,
    test_pred = pred_reg_ihe$test_pred$Nm,
    "ihe",
    var = "*Nm*"
  )

# COMBINE PLOTS  ---------------------------------------------------

layout <- "
AABBCCDD
EEFFGGHH
"

combo_plot <-
  roc_cat + conf_mat_cat + roc_ihe + conf_mat_ihe +
  op_cat_maxk + op_cat_nm + op_ihe_maxk + op_ihe_nm +
  plot_layout(design = layout) +
  plot_annotation(
    title = "**Model performance**",
    tag_levels = "A",
    tag_suffix = ")",
    theme = theme(plot.title = element_markdown(size = 20))
  )

ggsave(here("manuscript", "figures", "model_performance_fig.png"),
       combo_plot,
       width = 30,
       height = 15,
       units = "cm")

