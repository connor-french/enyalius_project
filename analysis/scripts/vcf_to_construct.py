import allel
import numpy as np

vcfs = snakemake.input

cs_out = snakemake.output


for i in range(len(vcfs)):
  callset = allel.read_vcf(vcfs[i])
  
  # convert to allele frequencies
  gt = allel.GenotypeArray(callset["calldata/GT"])
  gc = gt.count_alleles()
  gf = gc.to_frequencies()
  
  # only retain segregating alleles
  is_seg = gc.is_segregating()
  gf_seg = gf.compress(is_seg)
  
  # transpose so samples are rows and loci are columns
  gf_trans = np.transpose(gf)
  
  # write to file
  np.savetxt(cs_out[i], gf_trans)
