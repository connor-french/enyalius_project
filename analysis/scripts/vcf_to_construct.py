import allel
import numpy as np

vcfs = snakemake.input

cs_out = snakemake.output

for i in range(len(vcfs)):
  callset = allel.read_vcf(vcfs[i])

  # convert to allele counts and counts of the alternate allele
  gt = allel.GenotypeArray(callset["calldata/GT"])

  # get the allele counts matrix for all individuals
  gc = gt.count_alleles(max_allele = 1)

  # only retain segregating alleles
  is_seg = gc.is_segregating()
  gt_seg = gt.compress(is_seg, axis=0)
  gc_seg = gc.compress(is_seg)

  # alt allele counts matrix
  # needed for finding unlinked SNPs
  gn = gt_seg.to_n_alt(fill=-1)

  # filter for unlinked loci
  # locate unlinked SNPs
  loc_unlinked = allel.locate_unlinked(gn, threshold=0.1)

  # select unlinked SNPs
  gc_unlinked = gc_seg[loc_unlinked]
  gn_unlinked = gn[loc_unlinked]

  # remove singletons
  not_singleton = ~gc_unlinked.is_singleton(allele = 1)
  gn_unlinked = gn_unlinked[not_singleton]

  # convert to frequencies
  # na values will be -9
  gf = gn_unlinked / 2
  gf[gf < 0] = -9
  
  # transpose so samples are rows and loci are columns
  gf_trans = np.transpose(gf)
  
  # write to file
  np.savetxt(cs_out[i], gf_trans, fmt='%1.4f')
