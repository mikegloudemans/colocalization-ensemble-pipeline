import sys
import subprocess
import pandas as pd
from shutil import copyfile
import os
import math

(in_file,
	out_file) = sys.argv[1:]

trait1 = in_file.split("TRAIT1-")[1].split(".TRAIT2")[0]
trait2 = in_file.split("TRAIT2-")[1].split("/locus.")[0]

# Min cap for p-values so we're not taking a log of 0
min_allowable = 1e-150

def main():

	# Load data
	data = pd.read_csv(in_file, sep="\t")

	N_snps = data.shape[0]

	chrom = data["chr_trait1"].iloc[0]
	snp_pos = data["seed_pos"].iloc[0]
	
	baseline1, baseline2, baseline3, baseline4 = get_baseline_scores(data, chrom, snp_pos)

	min_pvalue_trait1 = min(list(data["pvalue_trait1"]))
	min_pvalue_trait2 = min(list(data["pvalue_trait2"]))
	feature1 = data["feature_trait1"].iloc[0]
	feature2 = data["feature_trait2"].iloc[0]
	
	# TODO: Write output
	with open(out_file, "a") as a:
		a.write(f"{chrom}\t{snp_pos}\t{trait1}\t{trait2}\t{min_pvalue_trait1}\t{min_pvalue_trait2}\t{feature1}\t{feature2}\t{N_snps}\t{baseline1}\t{baseline2}\t{baseline3}\t{baseline4}\n")

def get_baseline_scores(data, chrom, snp_pos):
	
	combined = data.copy()

	pval1 = 2
	pval2 = 2 # Will combine p-values from both traits (multiplied)
	pval3 = 2 # Will combine p-values from the GWAS and eQTL (max)
	pval4 = 2 # Will combine p-values from the GWAS and eQTL (log-ish)
	hit = combined[(combined["chr_trait1"] == chrom) & (combined["snp_pos"] == snp_pos)]
	if hit.shape[0] == 1:
		# If the locus is in the data, return its p-value
	
		trait2_pval = hit["pvalue_trait2"].iloc[0] + min_allowable
		trait1_pval = hit["pvalue_trait1"].iloc[0] + min_allowable

		pval1 = trait2_pval
		pval2 = trait2_pval * trait1_pval # Combine the two pvalues
		pval3 = max(trait2_pval, trait1_pval)	# The worse of the two values
		pval4 = (-1 / math.log10(trait2_pval)) / (-1 * math.log10(trait1_pval)) # Some more moderated weighting I guess... 

	else:
		# Otherwise return best p-value within 10000 bp
		region = combined[(combined["chr_trait1"] == chrom) & (abs(combined["snp_pos"] - snp_pos) <= 10000)]
		if len(list(region["pvalue_trait2"])) == 0:
			return((pval1, pval2, pval3, pval4))
			
		trait2_pval = min(region["pvalue_trait2"]) + min_allowable

		pval1 = trait2_pval
		print(region.head())
		idx = [i for i in range(len(list(region["pvalue_trait2"]))) if (math.log(list(region["pvalue_trait2"])[i]) - math.log(trait2_pval) - min_allowable) < 0.0000001][0]
		trait1_pval = list(region["pvalue_trait1"])[idx] + min_allowable

		pval2 = trait2_pval * trait1_pval 
		pval3 = max(trait1_pval, trait2_pval)
		pval4 = (-1 / math.log10(trait2_pval)) / (-1 * math.log10(trait1_pval))

	return((pval1, pval2, pval3, pval4))

if __name__ == "__main__":
	main()

