import allel
import numpy as np
import allel
import numpy as np
import sys
import argparse
import allel
import numpy as np
import sys

# Read the VCF file and specify the fields to extract

# Create an argument parser
parser = argparse.ArgumentParser(description='Calculate average depth from VCF file.')

# Add the vcf_path argument
parser.add_argument('vcf_path', help='Path to the VCF file')

# Parse the command-line arguments
args = parser.parse_args()

# Access the vcf_path argument value
# vcf_path = "assembly/full/clust92_outfiles/clust92.vcf"  # this is the path to the vcf used in the paper
vcf_path = args.vcf_path

# Read the VCF file and specify the fields to extract
callset = allel.read_vcf(vcf_path, fields=['calldata/DP'])
nonzero_vars = callset["calldata/DP"][np.ma.nonzero(callset["calldata/DP"])]
avg_depth = np.mean(nonzero_vars)
print(avg_depth)

