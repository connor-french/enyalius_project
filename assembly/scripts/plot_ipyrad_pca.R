# functions to plot PCA output


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


  # create a dataframe for plotting
  pc_df <- pca_obj %>%
    # filter for the desired species
    filter(name == species_name) %>%
    bind_cols(locs_df %>% filter(species == species_name))



  ggplot(pc_df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[color_var]])) +
    geom_point(alpha = 0.9, size = 3, shape = 21) +
    scale_fill_viridis_d() +
    labs(title = paste0("E. ", species_name)) +
    theme_bw()
}


# function to plot the cumulative variances of the first 8 PCs for each clustering threshold
# var_obj should be a data frame with a name column and var column
# species_name needs to be the species name
# clust_thresh needs to be a vector of cluster thresholds

plot_cumvar <- function(var_obj, species_name, clust_thresh) {

  vardf <- var_obj %>%
    filter(str_detect(name, species_name)) %>%
    mutate(clust_thresh = str_split_fixed(name, "_", n = 2)[,2]) %>%
    group_by(name, clust_thresh) %>%
    slice(1:8) %>%
    summarize(cumvar = sum(var)) %>%
    ungroup()

  p <- ggplot(vardf, aes(x = clust_thresh, y = cumvar)) +
    geom_point(size = 3) +
    labs(y = "Cumulative variance in first 8 PCs") +
    theme_bw()

  return(p)

}


