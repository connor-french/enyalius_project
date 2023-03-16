# script to branch assemblies to different clustering threshold values and run
import ipyrad as ip

# load the initial assembly
data = ip.load_json(snakemake.input.p1)

# the clustering threshold values
cb = snakemake.params.cb

for i in cb:
  b = data.branch('clust' + str(i)) # create branches
  b.set_params("clust_threshold", i * 0.01) # make new params files for each branch with the new param
  b.run("34567", auto=True, force=True)
