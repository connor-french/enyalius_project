# there are inland individuals that need to be removed prior to simulations for E. catenatus. 
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import rasterio
from rasterio.mask import mask

def crop_rast_buffer(rast, coords, buffer_dist, out_path):

  xy = [Point(xy) for xy in coords]
  gdf = gpd.GeoDataFrame(geometry=xy)
  # create convex hull
  convex_hull = gdf.unary_union.convex_hull

  # add buffer in degrees
  buffered_polygon = convex_hull.buffer(buffer_dist)

    # Load the raster file
  with rast as src:
      # Crop the raster with the buffered polygon
      out_image, out_transform = mask(src, [buffered_polygon], crop=True)
      out_meta = src.meta

  # Update the metadata to reflect the new image size
  out_meta.update({"driver": "GTiff",
                  "height": out_image.shape[1],
                  "width": out_image.shape[2],
                  "transform": out_transform})

  # Write the cropped image to a new file
  with rasterio.open(out_path, "w", **out_meta) as dest:
      dest.write(out_image)

filepath = "output/sdm_projections/projections_cat.tif"

r = rasterio.open(filepath)

localities = pd.read_csv("data/enyalius_locs_genetics.csv")

localities = localities[localities["species"] == "catenatus"]

# remove inland individuals
inland_inds = ["cat_jeq_MTR17108", "cat_jeq_MTR17172", "cat_jeq_MTR17377", "cat_mig_MTR19889", "cat_mig_RPD188", "cat_rem_MTR19915", "cat_rem_RPD128"]

localities = localities[~localities['id_code'].isin(inland_inds)]

# isolate the coordinates as a list of tuples
coords = list(zip(localities.longitude, localities.latitude))


crop_rast_buffer(r, coords=coords, buffer_dist = 2, out_path = "output/sdm_projections/projections_cat_inland-removed.tif")
