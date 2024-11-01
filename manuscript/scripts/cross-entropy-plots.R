# cross-entropy plots
library(tidyverse)
library(patchwork)
library(here)

d_ihe <- read_csv(here("manuscript", "data", "cross-entropy-scores-ihe.csv")) %>%
  na.omit()

plot_ce <- function(d, a = "a1") {
  d %>%
    filter(alpha == a) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(fill = "skyblue", color = "black", shape = 21, size = 3) +
    labs(x = "K", y = "Cross-entropy", title = a) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
}


# ihe plots ---------------------------------------------------------------

pa1_ihe <- plot_ce(d_ihe, "a1")
pa10_ihe <- plot_ce(d_ihe, "a10")
pa100_ihe <- plot_ce(d_ihe, "a100")

ce_ihe <- pa1_ihe /
  pa10_ihe /
  pa100_ihe

ggsave(here("manuscript", "figures", "cross-entropy-ihe.png"), ce_ihe, width = 10, height = 8)

# cat plots ---------------------------------------------------------------

d_cat <- read_csv(here("manuscript", "data", "cross-entropy-scores-cat.csv")) %>%
  na.omit()

pa1_cat <- plot_ce(d_cat, "a1")
pa10_cat <- plot_ce(d_cat, "a10")
pa100_cat <- plot_ce(d_cat, "a100")

ce_cat <- pa1_cat /
  pa10_cat /
  pa100_cat

ggsave(here("manuscript", "figures", "cross-entropy-cat.png"), ce_cat, width = 10, height = 8)





