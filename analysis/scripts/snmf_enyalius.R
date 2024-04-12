### This script is interactive.
### I'm reusing this for each species rather than making a new script for each species or making this into a function.

### load packages. LEA can be installed on bioconductor
library(LEA)
library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(glue)
library(here)

#POPsutilities contains helper functions to facilitate plotting and data manipulation. Downloaded from sNMF website
source(here("analysis", "scripts", "POPSutilities.R"))


#function to format data for plotting
get_plot_data <- function(q_matrix, samples, best_k){
  tbl <- as_tibble(q_matrix)

  colnames(tbl) <- c(paste0("pop_", 1:best_k))

  data <- tbl %>%
    mutate(id = samples) %>%
    gather('pop', 'prob', colnames(tbl)) %>%
    group_by(id) %>%
    mutate(likely_assignment = pop[which.max(prob)],
           assignment_prob = max(prob)) %>%
    arrange(likely_assignment, desc(assignment_prob)) %>%
    ungroup() %>%
    mutate(id = forcats::fct_inorder(factor(id)))

  return(data)
}

#plot sNMF function
plot_snmf <- function(plot_data) {
  snmf_plot <- ggplot(plot_data, aes(id, prob, fill = pop)) +
    geom_col() +
    scale_fill_viridis_d() +
    facet_grid(~likely_assignment, scales = 'free', space = 'free') +
    theme(axis.ticks = element_blank(), axis.text = element_blank(),
          axis.line = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), legend.position = "none",
          strip.text.x = element_text(size = 15, colour = "black", face = "bold"),
          strip.background = element_rect(fill = "transparent", colour = "black", linetype = "solid"),
          panel.background = element_rect(fill = "transparent"))

  return(snmf_plot)
}

# function to run snmf for different alpha values and write results to file
run_snmf <- function(species,  sample_names, out_dir){

  for (a in c(1, 10, 100)) {
    # copy the geno file so it has the alpha number at the end. sNMF creates its project folders from this filename, so it has to be unique across alpha values
    geno_path <- file.path(out_dir, glue("{species}.geno"))
    geno_copy_path <- file.path(out_dir, glue("{species}_a{a}.geno"))

    file.copy(geno_path, geno_copy_path)

    snmf_obj <-
      snmf(
        geno_copy_path,
        K = 1:5,
        project = "new",
        alpha = a,
        entropy = TRUE,
        repetitions = 20,
        CPU = 6,
        ploidy = 2
      )

    # plot the cross-entropy criterion
    ce_out <- file.path(out_dir, glue("cross-entropy_{species}_a{a}.png"))
    png(ce_out, height = 10, width = 15, units = "cm", res = 300, bg = "transparent")
    plot(snmf_obj, cex = 1.0, col = "lightblue", pch = 19, main = glue("a{a}"))
    dev.off()

    for (k in 1:5){
      # get the best run for the K value
      ce <- cross.entropy(snmf_obj, K = k)
      best_run <- which.min(ce)

      # get the Q matrix for plotting
      q_mat <- Q(snmf_obj, K = k, run = best_run)

      # plot the Q matrix and write both data and plot to file
      plot_data <- get_plot_data(q_mat, samples = sample_names, best_k = k)

      write_csv(plot_data, file.path(out_dir, glue("pop_assignment_a{a}_k{k}.csv")))

      snmf_plot <- plot_snmf(plot_data)

      snmf_out <- file.path(out_dir, glue("snmf-barplot_{species}_a{a}_k{k}.png"))

      ggsave(snmf_out, snmf_plot, width = 15, height = 5, units = "cm")
    }

    # remove the copied geno file
    file.remove(geno_copy_path)

    print(glue("Finished sNMF for a{a}!"))

  }

}


# e. perditus -------------------------------------------------------------

# get sample names
vcf_small_per <- read.vcfR(here("analysis", "data", "vcfs", "clust92_perditus_0.6_nosingletons.vcf"), nrows = 5)
sample_names_per <- colnames(vcf_small_per@gt)[-1]

run_snmf(
  species = "perditus",
  sample_names = sample_names_per,
  out_dir = here("analysis", "output", "snmf", "perditus")
)




