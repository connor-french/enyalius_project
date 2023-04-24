import ipyrad.analysis as ipa
import pandas as pd

def run_pca(species, sp_map, snps_path="", axes_out="", var_out="", is_clust = False,  mincov = 0.75, impute_method = "sample"): 
  
  if species == "all":
    imap_sp = sp_map
  else:
    imap_sp = {species: sp_map[species]}
  
  
  pca_sp = ipa.pca(
    data = snps_path,
    imap = imap_sp,
    minmap = minmap,
    mincov = mincov,
    impute_method = impute_method
  )
  
  pca_sp.run()

  # define column names
  col_names = []
  ncols = pca_sp.pcaxes[0].shape[0]
  for i in range(ncols):
    axnum = i + 1
    n = "PC_{0}".format(axnum)
    col_names.append(n)
  
  # convert pc axes array to a pandas dataframe and write out
  pca_df = pd.DataFrame(pca_sp.pcaxes[0], columns = col_names)
  
  pca_df.to_csv(axes_out)
  
  # if this is a PCA for cluster evaluation, convert variances to a dataframe and write out
  if is_clust:
    var_df = pd.DataFrame(pca_sp.variances[0], columns = ["var"])
    var_df.to_csv(var_out)
  
# create an imap for running PCAs in ipyrad
def create_imap(locs_df, splist):
  id_df = locs_df[["species", "id_code"]]
  id_df = id_df[id_df["species"].isin(splist)]
  
  imap = {}
  
  for i in splist:
    imap[i] = id_df["id_code"][id_df["species"] == i].tolist()
    
  return imap



# create an imap for conducting PCA per species
locs_nobad = pd.read_csv("processed_localities/locs_nobad.csv")

species_list = snakemake.params.splist

imap_sp = create_imap(locs_df=locs_nobad, splist=species_list[:-1])


# require that 50% of samples have data in each group
minmap = {i: 0.5 for i in imap_sp}


# run pca for batch effects. This is a single clustering threshold for each species and all species
for i in range(len(species_list)):
  run_pca(species_list[i], snps_path = snakemake.input[3], axes_out = snakemake.params.batch_out[i], is_clust = False, sp_map = imap_sp, mincov = 0.75, impute_method = "sample")

# run pca for cluster thresholds
pca_out = snakemake.params.pca_out
clust_thresh = snakemake.params.clust_thresh

for i in range(len(species_list)):
  for j in range(len(clust_thresh)):
    ao = pca_out[0] + "pcs_" + species_list[i] + "_" + str(clust_thresh[j]) + ".csv"
    vo = pca_out[0] + "var_" + species_list[i] + "_" + str(clust_thresh[j]) + ".csv"
    run_pca(species_list[i], snps_path = snakemake.input[j], axes_out = ao, var_out = vo, is_clust = True, sp_map = imap_sp, mincov = 0.75, impute_method = "sample")
