# function to convert sharing matrix to a vector for performing correlations
# share_obj is a tibble where the
# ingroup_names are the names of ingroup individual names
## they look like

share_to_vec <- function(share_obj, ingroup_names) {

  sdf <- share_obj %>%
    rename(sample = ...1) %>%
    filter(sample %in% ingroup_names) %>%
    select(sample,
           any_of(ingroup_names)) %>%
    # rearrange rows according to ingroup_names so the columns and rows are symmetrical
    arrange(factor(sample, levels = ingroup_names))

  sdf_mat <- as.matrix(sdf %>% select(-sample))

  sdf_vec <- sdf_mat[lower.tri(sdf_mat)]

  return(sdf_vec)
}
