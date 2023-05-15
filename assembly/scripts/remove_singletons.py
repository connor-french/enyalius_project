# short script to remove singletons from VCFs
import subprocess

in_vcfs = snakemake.input

out_vcfs = snakemake.output

for i in range(len(in_vcfs)): 
  ivcf = in_vcfs[i]
  ovcf = out_vcfs[i]
  # output singletons
  subprocess.run(f"vcftools --vcf {ivcf} --singletons", shell = True)
  # trim the singletons file to only include the first two columns without headers (it's what vcftools likes)
  subprocess.run("awk -F'\t' -v OFS='\t' 'NR>1{print $1,$2}' out.singletons > clean.singletons", shell = True)
  # exclude singleton positions
  subprocess.run(f"vcftools --vcf {ivcf} --exclude-positions clean.singletons --recode --recode-INFO-all --stdout > {ovcf}", shell = True)
