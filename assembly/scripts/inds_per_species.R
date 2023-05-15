# extract individuals per species and write individual names to file
library(readr)
library(dplyr)
library(purrr)

f <- read_csv(snakemake@input[["csv"]])

bad_inds <- unlist(snakemake@params[["bad_inds"]])

f2 <- filter(f, !(id_code %in% bad_inds))

# group by species and split
f2_split <- f2 %>%
  group_split(species)

names(f2_split) <- map_chr(f2_split, \(x) unique(x$species))

# select the id_code column
f2_ids <- f2_split %>%
  map(\(x) pull(x, id_code))

## write to text file without header
map2(f2_ids, names(f2_ids), \(x, y) write_lines(x, file = paste0(y, "_inds.txt")))
