library(readr)
library(dplyr)

locs = read_csv(snakemake@input[[1]])
bad_inds = snakemake@config[[1]]

locs_filt <- locs %>%
  filter(!(id_code %in% bad_inds))

write_csv(locs_filt, snakemake@output[[1]])
