#!/bin/bash

#SBATCH --job-name=catenatus_sims     ## Name of the job
#SBATCH --partition=production
#SBATCH --ntasks=40             ## Number of tasks (analyses) to run
#SBATCH --cpus-per-task=1      ## The number of threads the code will use
#SBATCH --mem-per-cpu=8000M     ## Real memory(MB) per CPU required by the job.
#SBATCH --output=catenatus_output.o
#SBATCH --error=catenatus_error.e
#SBATCH --mail-type=ALL,TIME_LIMIT_50
#SBATCH --mail-user=cfrench@gradcenter.cuny.edu

# set wd to my user folder
cd /global/u/connor_gc_12
# source my bashrc so I can use micromamba environment
source ~/.bashrc
# activate my enyalius-sims micromamba environment
micromamba activate enyalius-sims
# run the script
python3 demography_simulations_cat.py

