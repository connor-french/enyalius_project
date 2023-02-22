# function to read in the _stats.txt files in the _outfiles/ folder for and ipyrad assembly

read_ipyrad_stats <- function(path) {
  t <- read_lines(path)
  ss <- t[grep("## Final Sample stats summary", t):length(t)] # extract the stats summary section
  d <- grep("## Alignment matrix statistics:", ss) - 3 # this is the end of the section, which we don't want
  ss <- ss[1:d]

  # convert the strings into a tibble
  s_full <- read_table(
    ss,
    skip = 2,
    col_names = c(
      "id_code",
      "state",
      "reads_raw",
      "reads_passed_filter",
      "clusters_total",
      "clusters_hidepth",
      "hetero_est",
      "error_est",
      "reads_consens",
      "loci_in_assembly"
    )
  ) %>%
    mutate(species = str_extract(id_code, "^[^_]+(?=_)"))
}

