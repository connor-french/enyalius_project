# function to obtain a matrix of shared alleles among individuals to visualize patterns of missingness (i.e. less sharing = more missing)
from ipyrad.analysis.sharing import Sharing
import pandas

def run_sharing(data, outpath):
  share_cl = Sharing(data = data)
  share_cl.run()
  share_cl.sharing_matrix.to_csv(outpath)
  print("Written to {0}".format(outpath))

for i in range(len(snakemake.input)):
  run_sharing(snakemake.input[i], snakemake.output[i])
