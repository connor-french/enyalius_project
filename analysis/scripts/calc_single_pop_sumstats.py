import allel
import numpy as np
import scipy.spatial.distance as spd

vcfs = snakemake.input
taj_d_out = snakemake.output["taj_d"]
site_pi_out = snakemake.output["site_pi"]
ind_pi_out = snakemake.output["ind_pi"]
wat_theta_out = snakemake.output["wat_theta"]
sfs_out = snakemake.output["sfs"]
num_sites_out = snakemake.output["num_sites"]

for i in range(len(vcfs)):
  callset = allel.read_vcf(vcfs[i])
  
  # convert to allele counts and counts of the alternate allele
  gt = allel.GenotypeArray(callset["calldata/GT"])
  gc = gt.count_alleles()
  
  # only retain segregating alleles
  is_seg = gc.is_segregating()
  gt_seg = gt.compress(is_seg, axis=0)
  gc_seg = gc.compress(is_seg)
  gn = gt_seg.to_n_alt()

  # filter for unlinked loci
  # locate unlinked SNPs
  loc_unlinked = allel.locate_unlinked(gn, threshold=0.1)

  # select unlinked SNPs
  gc_unlinked = gc_seg[loc_unlinked]
  gn_unlinked = gn[loc_unlinked]
  
  # Tajima's D in moving windows of 100 bp (non-overlapping)
  taj_d = allel.moving_tajima_d(gc_unlinked, size=100)
  
  # per-site pi (mean pairwise nucleotide differences across SNPs)
  site_pi = allel.mean_pairwise_difference(gc_unlinked)
  
  # per-ind pi (number of nucleotide differences across individuals)
  ind_pi = allel.pairwise_distance(gn_unlinked, metric="euclidean")
  ## convert to a square matrix
  ind_pi_mat = spd.squareform(ind_pi)
  
  # watterson's theta in moving windows of 100 bp (non-overlapping)
  # don't need unlinked loci and you need non-segregating sites
  end = gc.shape[0]
  pos = np.arange(0, end, 1)
  
  wat_theta = allel.windowed_watterson_theta(pos, gc, size=100)
  
  # folded site frequency spectrum
  sfs = allel.sfs_folded(gc_unlinked)
  
  # number of sites
  num_sites = gn_unlinked.shape[0]
  
  # write to file
  
  ## site-pi
  np.savetxt(site_pi_out[i], site_pi)
  
  ## ind-pi
  np.savetxt(ind_pi_out[i], ind_pi_mat)
  
  ## Watterson's theta
  np.savetxt(wat_theta_out[i], wat_theta[0])
  
  ## tajima's d
  np.savetxt(taj_d_out[i], taj_d)
  
  ## sfs
  np.savetxt(sfs_out[i], sfs)
  
  ## number of sites
  with open(num_sites_out[i], 'w') as f:
    f.write(str(num_sites))

