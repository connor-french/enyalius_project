README
================

# Enyalius project

Code and most data associated with the *Enyalius* comparative
phylogeography project. The README right now is just a file guide.

Directory and file guide:

-   `analysis/`- where I house all analysis-related files. It’s
    organized into subfolders:

    -   `data/`- mostly raw data files used for analysis

    -   `atlantic_forest/atlantic_forest.geojson`- shapefile of the
        Atlantic Forest boundaries. I should figure out the provenance
        of this file- I don’t remember where I got it.

    -   `brazil_clim/`- bioclimate variables (`.tif` rasters) refined
        for Brazil, from [Ramoni-Perazzi et
        al. 2022](https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/joc.7325).
        Downloaded on 2022-02-25.

    -   `chelsa_trace21k/`- CHELSA bioclimatic variables from the LGM

    -   `forest_cover/`- forest cover rasters (classes 1-4) from [Tuanmu
        and Jetz
        2014](https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.12182).
        The data was housed at
        [earthenv.org/landcover](https://www.earthenv.org/landcover).
        Downloaded on 2022-02-26.

    -   `enyalius_locs.csv`- table of enyalius localities and metadata

    -   `Laura_ALL_ms_FINAL_80complete.snps.phy`- SNPs file I’m using
        until I finalize my own assemblies of the data

    -   `output/`- output files from analyses. This folder is labile, so
        will change regularly.

        -   `cropped_predictors/`- cropped bioclims and forest cover
            variables for SDMs

        -   `sdm_models/`- SDM models for each *Enyalius* species

        -   `sdm_projections/`- SDM projections for each species

        -   `sdm_response_curves/`- response curves for environmental
            variables used in the SDMs

        -   `thinned_localities/`- localites for each species after
            spatial thinning. Saved as geojson files

    -   `reports/`- Rmarkdown and Quarto reports for each analysis

        -   `assembly_exploration.Rmd`- interactive document to explore
            assembly statistics. Runs on shiny.

        -   `exploratory_phylogenetic_analysis.Rmd`- some quick trees to
            investigate if any individuals don’t map to their species or
            are otherwise wonky

        -   `explore-slim-output.qmd`- exploring SLiM output of spatial
            simulations that I’m toying around with. I used them for the
            2022 Evolution meeting

        -   `slim_map_small.qmd`- short script to aggregate and convert
            SDM rasters to SLiM’s format. They are at a lower resolution
            to facilitate speed

        -   `species_distribution_modeling.Rmd`- species distribution
            models of each species with some post-processing

        -   `spatial_thin_log.txt`- log file produced when running
            spThin to spatially thin localities

        -   `recap/`- msprime recapitation output that gets generated
            when processing SLiM output

        -   `*_files/`- static files that get generated when rendering
            `.qmd` reports

    -   `slim/`- experimenting with running spatial simulations in SLiM.
        These are modified from the Peter Ralph’s
        [nebria](https://github.com/petrelharp/nebria) github repo.

-   `assembly/`- files related to *Enyalius* RADseq assembly. This is
    going to be modified soon \* `man/`- folder to house documentation
    of functions if I turn this into a package
    
    -   `processed_localities/`- spreadsheets containing data about samples with locality data that were processed from the original raw CSV. 
        -   `locs_nobad.csv`- all individuals that had poor sequencing or otherwise had something wrong with them removed

-   `R/`- R scripts used for analysis

    -   `make_maps_fns.R`- functions to convert a raster to a SLiM map.
        Taken from Peter Ralph’s `nebria` repository
    -   `adhoc/` - functions used for a specific/single use that aren’t
        needed to reproduce results
        -   `conda_setup.R` - script I used to setup a conda environment
            for python package management. Users don’t need it because
            `renv` tracks python packages after initializing the conda
            repo. They just need to install miniconda or anaconda

        -   `evolution_2022_figures.R`- code I used to generate figures
            for the 2022 Evolution Conference in Cleveland

-   `renv/`- folder to house R and Python package info. Don’t touch
    anything here.
