# README


<figure>
<a href="https://doi.org/10.5281/zenodo.14029198"><img
src="https://zenodo.org/badge/460098004.svg" alt="DOI" /></a>
<figcaption>DOI</figcaption>
</figure>

Code and access to data associated with the *Enyalius* comparative
phylogeography project.

First, clone this repository.

To install python packages necessary for analysis, first install
[miniconda](https://docs.conda.io/en/latest/miniconda.html) if you do
not already have anaconda/miniconda. Then, create a conda environment
with the command:

`conda env create -f environment.yml`

To activate the environment, run:

`conda activate enyalius`

### Running assembly and analyses

1.  In the terminal, make sure you are in the main `enyalius/` directory
2.  Activate the conda environment with `conda activate enyalius`
3.  run `snakemake`. You may want to perform a dry run first with
    `snakemake -n`. See the [snakemake
    documentation](https://snakemake.readthedocs.io/en/v5.6.0/index.html)
    for more on using snakemake
4.  navigate to the `assembly/` directory
5.  run `snakemake`. CAUTION- this is very compute intensive. Make sure
    you have appropriate resources. You may want to run `snakemake` in
    parallel with the `snakemake -j` flag
6.  navigate to the `analysis/` folder
7.  Follow the steps in `species_distribution_modeling.md` to run the
    SDM analysis (may want to view the md file on Github)
8.  run `snakemake`
9.  Follow the steps in `empirical_sumstats.md` to calculate empirical
    summary statistics
10. run each file in the `scripts/` folder, following instructions in
    the file

### Directory and file guide

-   `analysis/`- where I house all analysis-related files. It’s
    organized into subfolders:

    -   `data/`- data files used for analysis

        -   `data/atlantic_forest/atlantic_forest.geojson`- shapefile of
            the Atlantic Forest boundaries

        -   `data/current_climate_chelsa/`- CHELSA bioclimate data for
            present-day

        -   `data/vcfs/`- filtered VCFs (output from assembly) for
            missing data exploration and final analysis

        -   `enyalius_locs.csv`- table of enyalius localities and
            metadata

        -   `*_inds.txt`- lists of individuals per species

    -   `empirical_sumstats.*` - notebook outlining the calculation of
        genetic summary statistics for the empirical samples

    -   `scripts/` - Python and R scripts used to conduct the analysis

    -   `single_pop_sumstats.*` - notebook exploring per-species genetic
        summary statistics calculated for different levels of missing
        data

    -   `species_distribution_modeling`.\* - notebook to run SDMs

    -   `Snakefile` - Snakemake config file for running all analyses

    -   `renv/` - R package management

    -   `renv.lock` - R package management

-   `assembly/`- files related to *Enyalius* RADseq assembly

    -   `assembly-qc*` - report outlining assembly quality control
    -   `config.yaml` - configuring file for Snakemake run
    -   `fastq/` - fastq files (are summarized in `multiqc.html`)
    -   `processed_localities/` - localities used for assembly filtering
    -   `scripts/` - scripts used for assembly
    -   `*inds.txt`

-   `R/`- R scripts used for analysis

    -   `make_maps_fns.R`- functions to convert a raster to a SLiM map.
        Taken from Peter Ralph’s `nebria` repository

-   `renv/`- folder to house R and Python package info. Don’t touch
    anything here.

-   `Snakefile` - snakemake config to move files around for analyses and
    assembly
