# for whatever reason, for loops don't want to return the output of run() to the loop, so I need to put it in a function
from ipyrad.analysis.sharing import Sharing
import pandas

def run_sharing(data, outpath):
  share_cl = Sharing(data = data)
  share_cl.run()
  share_cl.to_csv(outpath)
  print("Written to {0}".format(outpath))

for i in range(len(snakemake.input)):
  run_sharing(snakemake.input[i], snakemake.output[i])
