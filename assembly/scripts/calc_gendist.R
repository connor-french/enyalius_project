# calculate genetic distance in PC-space
calc_gendist <- function(pca_obj, num_pcs = 20) {

  # get the pc axes
  pcs <- pca_obj$pcaxes$`0`[,1:num_pcs]

  # convert pc axes dataframe to a matrix
  pc_dist <- as.matrix(dist(pcs, method = "euclidean"))

  # get the lower triangle of the matrix as a vector
  pc_dist_vec <- pc_dist[lower.tri(pc_dist)]

  return(pc_dist_vec)
}
