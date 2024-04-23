# sequencing summary
library(tidyverse)
library(here)

# final set of individuals after removing outliers
locs <- read_csv(here("analysis", "data", "enyalius_locs_genetics.csv")) %>%
  filter(!str_detect(id_code, "pic_"))


# assembly statistics
seqstats <- read_table(
  here("assembly", "full", "clust92_outfiles", "clust92_stats.txt"),
  col_names = c("id", "state",  "reads_raw",  "reads_passed_filter",  "clusters_total",  "clusters_hidepth",  "hetero_est",  "error_est",  "reads_consens",  "loci_in_assembly"),
  skip = 386,
  n_max = 161
) %>%
  filter(
    !str_detect(id, "pic_")
  )

# average number of reads pre-filter
seqstats %>%
  summarize(mean_reads = mean(reads_raw))

# average number of reads post-filter
seqstats %>%
  summarize(mean_reads = mean(reads_passed_filter))


# average of all stats
seqstats %>%
  summarize(across(where(is.numeric), mean)) %>% View()

