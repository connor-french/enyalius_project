# script to plot predictive intervals of empirical data and barplots of classification probability
library(tidyverse)
library(ggtext)
library(ggrepel)
library(glue)
library(here)

theme_set(theme_bw())
theme_update(plot.title = element_markdown(size = 16, face = "bold"),
             axis.title = element_text(size = 14),
             axis.text = element_text(size = 12))

# 0. functions ------------------------------------------------------------

read_emp_pred <- function(species, param, prior_range) {
  folder <- here("analysis", "output", "parameter_estimation", glue("reg_tuning_{species}"))
  df <- read_csv(here(folder, glue("emp_pred_{param}.csv"))) %>%
    mutate(species = species,
           prior_lower = prior_range[1],
           prior_upper = prior_range[2])

  return(df)
}

plot_pred_int <- function(df, param) {

  if (param == "max_k") {
    param_lab <- "Max deme size"
    ynudge <- 300
    rnd <- 0
  } else if (param == "Nm") {
    param_lab <- param
    ynudge <- 0.1
    rnd <- 2
  } else {
    param_lab <- param
    ynudge <- 1000
    rnd <- 0
  }

  df <- df %>%
    mutate(sp_long = case_when(
      species == "cat" ~ "E. catenatus",
      species == "ihe" ~ "E. iheringii",
      species == "per" ~ "E. perditus",
      .default = species
    ))

  p <- ggplot(df) +
    geom_segment(aes(
      x = sp_long,
      xend = sp_long,
      y = lower,
      yend = upper,
      color = species
    ),
    linewidth = 3) +
    geom_point(aes(x = sp_long, y = center, color = sp_long), size = 8) +
    scale_color_manual(values = c("#DE8E07", "#DE8E07", "#0A5DA2", "#0A5DA2")) +
    geom_hline(yintercept = df$prior_lower, linetype = "dashed") +
    geom_hline(yintercept = df$prior_upper, linetype = "dashed") +
    geom_label_repel(
      aes(
        x = sp_long,
        y = center,
        label = round(center, rnd)
      ),
      size = 5,
      nudge_x = 0.3,
      nudge_y = ynudge,
      segment.linetype = "dashed",
      segment.color = "darkgray",
      point.padding = 1.1
    ) +

    labs(y = param_lab, title = "95% prediction intervals") +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_markdown(face = "italic")
    )
  return(p)
}

plot_barplot <- function(species, prob_thresh) {

  prob_lin <- 1 - prob_thresh
  if (species == "cat") {
    sp_long <- "E. catenatus"
    fill_color <- "#DE8E07"
  } else if (species == "ihe") {
    sp_long <- "E. iheringii"
    fill_color <- "#0A5DA2"
  } else {
    sp_long <- "E. perditus"
    fill_color <- "darkgreen"
  }

  df <- tibble(
    species = species,
    prob = c(prob_thresh, prob_lin),
    label = c("Threshold", "Linear")
  )

  p <- ggplot(data = df) +
    geom_bar(aes(x = label, y = prob), stat = "identity", fill = fill_color, width = 0.6) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent) +
    labs(title = glue("Empirical classification probabilities for<br>*{sp_long}*</br>"), y = "Probability") +
    theme(axis.title.y = element_blank())

  return(p)
}

# 1. load data ------------------------------------------------------------

emp_pred_max_k <- map_df(c("cat", "ihe"), ~read_emp_pred(.x, "max_k", c(100, 5000)))
emp_pred_nm <- map_df(c("cat", "ihe"), ~read_emp_pred(.x, "Nm", c(0.1, 5)))

# 2. plot data ------------------------------------------------------------
p_max_k <- plot_pred_int(emp_pred_max_k, "max_k")
p_nm <- plot_pred_int(emp_pred_nm, "Nm") + theme(axis.title.y = element_markdown(face = "italic"))

p_cat <- plot_barplot("cat", 0.6)
p_ihe <- plot_barplot("ihe", 0.99)


(p_cat + p_ihe) / (free(p_max_k) + free(p_nm)) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ")"
  ) & theme(plot.tag = element_text(size = 16))

# 3. save plot ------------------------------------------------------------

ggsave(
  here("manuscript", "figures", "pred_int_fig.png"),
  width = 30,
  height = 20,
  units = "cm",
  dpi = 600
)








