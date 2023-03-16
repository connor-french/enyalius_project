def run_pca(species, snps_path, imap, mincov = 0.75, impute_method = "sample"): 
  
  if species == "all":
    imap_sp = imap
  else:
    imap_sp = {species: imap[species]}
  
  
  pca_sp = ipa.pca(
    data = snps_path,
    imap = imap_sp,
    minmap = minmap,
    mincov = mincov,
    impute_method = impute_method
  )
  
  pca_sp.run()
  
  return pca_sp
