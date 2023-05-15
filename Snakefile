configfile:
  "config.yaml"

rule all:
  input:
    expand("analysis/data/vcfs/clust92_{species}_{minspcov}_nosingletons.vcf", species = config["INLONG"], minspcov = config["MINSPCOV"]),
    "analysis/data/enyalius_locs_genetics.csv",
    expand("analysis/data/{species}_inds.txt", species = config["INLONG"])

# move singleton vcfs to analysis folder

rule move_singletons:
  input:
    expand("assembly/full/clust92_outfiles/clust92_{species}_{minspcov}_nosingletons.vcf", species = config["INLONG"], minspcov = config["MINSPCOV"])
  output:
    expand("analysis/data/vcfs/clust92_{species}_{minspcov}_nosingletons.vcf", species = config["INLONG"], minspcov = config["MINSPCOV"])
  shell:
    "mv assembly/full/clust92_outfiles/*nosingletons.vcf analysis/data/vcfs/"

# make sure I have the localities for genetics individuals who don't have sequencing problems for analysis
rule move_localities:
  input:
    "assembly/processed_localities/locs_nobad_final.csv"
  output:
    "analysis/data/enyalius_locs_genetics.csv"
  shell:
    "cp {input} {output}"

# need individual lists in order for genetic data visualizations

rule move_ind_lists:
  input:
    expand("assembly/{species}_inds.txt", species = config["INLONG"])
  output:
    expand("analysis/data/{species}_inds.txt", species = config["INLONG"])
  shell:
    "cp assembly/*_inds.txt analysis/data/"
