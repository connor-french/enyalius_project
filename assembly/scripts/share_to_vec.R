# function to convert sharing matrix to a vector for performing correlations
# share_obj is a python object output from the function Sharing in ipyrad.analysis
# pca_names are the names of individuals from the PCA performed in ipyrad.analysis
## they look like

share_to_vec <- function(share_obj, pca_names) {

  sdf <- share_obj$sharing_matrix %>%
    as_tibble(rownames = "sample") %>%
    filter(sample %in% pca_names) %>%
    select(sample,
           any_of(pca_names))

  sdf_mat <- as.matrix(sdf %>% select(-sample))

  sdf_vec <- sdf_mat[lower.tri(sdf_mat)]

  return(sdf_vec)
}
