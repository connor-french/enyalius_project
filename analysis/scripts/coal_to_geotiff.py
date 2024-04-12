import numpy as np
import os
import rasterio
import argparse
import warnings

'''
This script reads in per-deme diversity CSVs output by demographic_simulations_posteriors.py,
takes the per-deme average, and converts it to a geotiff for visualization in R.
'''

# Create an argument parser
parser = argparse.ArgumentParser(description='Convert per-deme diversity CSVs to geotiff.')

# Add arguments for folder path and reference raster path
parser.add_argument('--folder', type=str, help='Path to the folder containing the CSV files.')
parser.add_argument('--refraster', type=str, help='Path to the reference raster.')
parser.add_argument('--outraster', type=str, help='Path to the output raster.')

# Parse the command line arguments
args = parser.parse_args()

# Get the folder path and reference raster path from the arguments
folder_path = args.folder
reference_raster_path = args.refraster
out_raster_path = args.outraster

# Rest of the code...

# Get a list of all CSV files in the folder
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]

# Initialize an empty list to store the arrays
arrays = []

# Read each CSV file and convert it to a numpy array
for file in csv_files:
  file_path = os.path.join(folder_path, file)
  array = np.loadtxt(file_path, delimiter=',')
  arrays.append(array)

# Convert the list of arrays to a single ndarray
combined_array = np.stack(arrays)

# Replace -1 with NaN in the combined_array
combined_array[combined_array == -1] = np.nan

# Take the element-wise mean of the combined_array
# I expect RuntimeWarnings for empty mean slices, so I'm ignoring them
with warnings.catch_warnings():
  warnings.simplefilter("ignore", category=RuntimeWarning)
  mean_array = np.nanmean(combined_array, axis=0)

# Open the reference raster
with rasterio.open(reference_raster_path) as src:
  # Get the metadata from the reference raster
  meta = src.meta

# update meta to tell it what the nodata is and how many layers there are
meta["nodata"] = np.nan
meta["count"] = 1

# Create a new raster file to write the mean_array
with rasterio.open(out_raster_path, 'w', **meta) as dst:
  # Write the mean_array to the new raster file
  dst.write(mean_array, 1)
