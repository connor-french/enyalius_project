# function to plot the first 4 PCs for a PCA run by ipyrad
# pca_obj needs to be output from ipyrad.analysis.pca
# x and y are the plot axes and need to be the form "PC_1"
# locs_df is a data frame that contains the id_code column and any other metadata
# color_var is the variable that you want to color the points by. Has to be a variable from locs_df

plot_pca <- function(pca_obj,
                     x = "PC_1",
                     y = "PC_2",
                     locs_df = locs,
                     species_name,
                     color_var) {

  pc_df <- suppressMessages(bind_cols(pca_obj$pcaxes[["0"]][,1:4],
                                      pca_obj$names))

  colnames(pc_df) <- c(paste0("PC_", 1:4), "id_code")

  pc_df_cat <- left_join(
    pc_df,
    locs_df,
    by = "id_code"
  )

  ggplot(pc_df_cat, aes(x = .data[[x]], y = .data[[y]], fill = .data[[color_var]])) +
    geom_point(alpha = 0.9, size = 3, shape = 21) +
    scale_fill_viridis_d() +
    labs(title = paste0("E. ", species_name)) +
    theme_bw()
}
