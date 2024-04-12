# command line arguments used to run the script and convert the coal output to geotiff

## iheringii
python3 coal_to_geotiff.py \
--folder "../output/posterior_sims/post_ihe/post_csv" \
--refraster "../output/sdm_projections/projections_ihe.tif" \
--outraster "../output/gendiv_rasters/gendiv_ihe.tif" 

## catenatus
python3 coal_to_geotiff.py \
--folder "../output/posterior_sims/post_cat/post_csv" \
--refraster "../output/sdm_projections/projections_cat_inland-removed.tif" \
--outraster "../output/gendiv_rasters/gendiv_cat.tif" 