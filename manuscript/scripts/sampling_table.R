# script to create a sampling table for the manuscript
library(geobr)
library(readxl)
library(gt)
library(sf)
library(tidyverse)
library(here)
library(glue)


# 1. read in spreadsheet --------------------------------------------------
sdf <- read_excel(here("manuscript", "tables", "spreadsheets", "genetic_sampling_table.xlsx"))

# have to use a version without NAs in lat/longs to get municipality information. This is just the outgroups
sdf_noout <- sdf %>%
  filter(latitude != "NA") %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )


# 2. extract municipality information from coordinates --------------------

# download municipalities
mun <- read_municipality(code_muni = "all", year = 2020)

# extract municipality information
sdf_mun <-
  sdf_noout %>%
  # transform to same crs as municipalities
  st_transform(st_crs(mun)) %>%
  st_join(mun, join = st_within) %>%
  # transform back to 4326
  st_transform(4326) %>%
  # add in the removed outgroups
  bind_rows(sdf %>% filter(latitude == "NA") %>% mutate(latitude = NA_character_), .id = "outgroup") %>%
  mutate(
    outgroup = if_else(outgroup == "1", "Ingroup", "Outgroup"),
    # add in municipality and state information for outgroups
    name_muni = if_else(outgroup == "Outgroup", glue("**{municipality}**"), name_muni),
    name_state = if_else(outgroup == "Outgroup", glue("**{state}**"), name_state),
    species = if_else(outgroup == "Outgroup", glue("***E. {species}***"), glue("*E. {species}*")),
    sample_id = if_else(outgroup == "Outgroup", glue("**{sample_id}**"), sample_id)
    )


# 4. create supplementary table ---------------------------------------------------------
samp_table_supp <- sdf_mun %>%
  as_tibble() %>%
  select(
    "Sample ID" = sample_id,
    "Species" = species,
    "Latitude" = latitude,
    "Longitude" = longitude,
    "Municipality" = name_muni,
    "State" = name_state,
    outgroup
  ) %>%
  mutate(
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude),
  ) %>%
  # remove perditus
  filter(Species != "*E. perditus*") %>%
  # arrange by species, then municipality
  arrange(outgroup, Species, State, Municipality) %>%
  select(-outgroup) %>%
  gt() %>%
  fmt_number(
    columns = c(Latitude, Longitude),
    decimals = 2
  ) %>%
  fmt_markdown(columns = c(Species, `Sample ID`, Municipality, State)) %>%
  # add a footnote for latitude and longitude
  tab_footnote(
    footnote = "Latitude and longitude are in decimal degrees.",
    locations = cells_column_labels(columns = c(Latitude, Longitude))
  ) %>%
  cols_align(
    align = "center",
    columns = c(Latitude, Longitude, Municipality, State)
  ) %>%
  cols_width(
    c(Latitude, Longitude, Municipality, State) ~ px(100)
  )

gtsave(samp_table_supp, here("manuscript", "tables", "supp_sampling_table.tex"))
gtsave(samp_table_supp, here("manuscript", "tables", "supp_sampling_table.docx"))

