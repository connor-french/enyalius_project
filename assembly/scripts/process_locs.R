library(here)
library(readr)
library(dplyr)

low_cov_inds <- c("cat_jibo_MTR18418", "cat_jeq_MTR17278", "cat_uba_PEU545", "cat_jib_LG1377", "cat_itc_MTR11831", "cat_itc_MTR11833", "cat_ruib_MTR18434", "cat_jibo_LG1377", "cat_jeq_MTR17431", "ihe_saop_CTMZ13066", "ihe_ipor_CTMZ11224", "ihe_nhae_CTMZ13229", "ihe_smig_CTMZ03618", "ihe_nhae_CTMZ13228", "per_petrop_MTR22755", "per_saop_UNIBAN0838", "per_petr_MTR22755", "per_saop_UNIBAN0883", "pic_tra_MTR13499", "pic_tra_MTR13500")

# individuals found in the missingness analysis to contain low sequencing coverage and low numbers of loci
bad_inds <- c("per_bana_CTMZ03943", "cat_cam_PEU322")

# read in localities
locs <- read_csv(here("analysis", "data", "enyalius_locs.csv")) %>%
  ## filter out low coverage individuals
  filter(!(id_code %in% low_cov_inds),
         is.na(obs) | obs != "Excluded low coverage")

# remove bad inds for analyses post-missingness assessment
locs_nobad <- locs %>%
  filter(!(id_code %in% bad_inds))

# write to file
write_csv(locs_nobad, here("assembly", "processed_localities", "locs_nobad.csv"))
