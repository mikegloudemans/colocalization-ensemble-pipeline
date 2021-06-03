# Author: Mike Gloudemans
#
# preprocess.py
#
# Tools for loading summary statistics.
#

import subprocess
import sys
import gzip
import pandas as pd
from io import StringIO
from scipy import stats

out_file = sys.argv[1]
window = int(sys.argv[2])
source_file = sys.argv[3]
lookup_file = sys.argv[4]
chrom = sys.argv[5]
pos = int(sys.argv[6])
source_trait = sys.argv[7]
lookup_trait = sys.argv[8]

def main():

	source_data = get_sumstats(source_file, chrom, pos, "_source", source_trait)
	lookup_data = get_sumstats(lookup_file, chrom, pos, "_lookup", lookup_trait)

	if isinstance(source_data, str) or source_data.shape[0] == 0:
		print(f"No overlapping data for {source_file} {lookup_file} {chrom} {pos} {source_trait} {lookup_trait}")
		return
	if isinstance(lookup_data, str) or lookup_data.shape[0] == 0:
		print(f"No overlapping data for {source_file} {lookup_file} {chrom} {pos} {source_trait} {lookup_trait}")
		return

	print(source_data.head(20))
	print(lookup_data.head(20))
	merge_data = combine_sumstats(source_data, lookup_data)
	if isinstance(merge_data, str) or merge_data.shape[0] == 0:
		print(f"No overlapping data for {source_file} {lookup_file} {chrom} {pos} {source_trait} {lookup_trait}")
		return

	merge_data.to_csv(out_file, sep="\t", index=False, header=True)

def get_sumstats(trait_file, chrom, pos, suffix, trait="none"):

	# Get data using tabix

	with gzip.open(trait_file, 'rb') as f:
		header = f.readline().decode('utf-8')

	min_pos = pos-window
	max_pos = pos+window

	data = subprocess.run(f"tabix {trait_file} chr{chrom}:{min_pos}-{max_pos}".split(), capture_output=True).stdout.decode('utf-8') + \
		subprocess.run(f"tabix {trait_file} {chrom}:{min_pos}-{max_pos}".split(), capture_output=True).stdout.decode('utf-8')

	table = pd.read_csv(StringIO(header + data), sep="\t")

	# Filter trait if needed
	if trait != "none":
		if "trait" in list(table.columns.values):
			table = table[table['trait'] == trait].copy()
		elif "feature" in list(table.columns.values):
			table = table[table['feature'] == trait].copy()
		elif "gene" in list(table.columns.values):
			table = table[table['gene'] == trait].copy()
		else:
			return("No trait column specified; trait filtering is impossible")
		
	if table.shape[0] == 0:
		return "No summary statistics found at this locus."

	table['pvalue'] = table['pvalue'].astype(float)

	# TODO: handle this exception

	if "ref" in list(table.columns.values):
		table['non_effect_allele'] = table['ref']		
		table['effect_allele'] = table['alt']		

	if "effect_allele" in list(table.columns.values):
		table['effect_allele'] = table["effect_allele"].str.upper()
		table['non_effect_allele'] = table["non_effect_allele"].str.upper()

	if 'zscore' in table.columns.values:
		table = table.rename(index=str, columns={"zscore": "ZSCORE"})

	# Derive zscore from beta and se if needed, or from p-value and direction
	if 'beta' in table and 'se' in table:
		table['ZSCORE'] = table['beta'] / table['se']
		table = table.fillna({'beta': 0})
		table['ZSCORE'] = table['ZSCORE'].fillna(0)

	else: 
		table = table[~table["pvalue"].isna()]

		# Need to cap it at z-score of 40 for outrageous p-values (like with AMD / RPE stuff)
		if "direction" in table:
			table['ZSCORE'] = pd.Series([min(x, 40) for x in stats.norm.isf(table["pvalue"] / 2)], index=table.index) * (2*(table["direction"] == "+")-1)

	# Now derive beta and se from zscore if needed
	if 'ZSCORE' in table.columns.values:
		if 'se' not in table.columns.values or 'beta' not in table.columns.values:
			table['beta'] = table['ZSCORE']
			table['se'] = 1

	# Tag with suffix for when the merge happens	
	if "beta" in table:
		table = table.rename(index=str, columns={"beta": "beta" + suffix})
	if "se" in table:
		table = table.rename(index=str, columns={"se": "se" + suffix})

	table['snp_pos'] = table['snp_pos'].astype(int)
	
	if min(table['pvalue']) < 1e-150:
		table['pvalue'] = table['pvalue'].apply(lambda x: max(x, 1e-150))

	# Sometimes the GWAS SNP is outside of the range of eQTLs tested for a certain
	# gene, or on the outside fringe of the range. If this is the case, then skip it.
	if pos > max(table['snp_pos']) - 50000 or pos < min(table['snp_pos']) + 50000:
		return "Sumstats interval too small."
	
	return table

def combine_sumstats(source_data, lookup_data):

	combined = pd.merge(source_data, lookup_data, on="snp_pos", suffixes=("_source", "_lookup"))
   
	# Remove all positions that appear multiple times in the GWAS table.
	dup_counts = {}
	for pos in combined['snp_pos']:
			dup_counts[pos] = dup_counts.get(pos, 0) + 1

	combined['dup_counts'] = [dup_counts[pos] for pos in combined['snp_pos']]
	combined = combined[combined['dup_counts'] == 1]

	# Check to make sure there are SNPs remaining; if not, just move on
	# to next gene.
	if combined.shape[0] == 0:
		return "No overlapping SNPs in eQTL and GWAS"

	# Check to make sure we still have significant GWAS hits and eQTLs, if desired
	#if "screening_thresholds" in settings and "gwas" in settings["screening_thresholds"]:
	#	if min(combined['pvalue_gwas']) > settings["screening_thresholds"]["gwas"]:
	#		return "No significant GWAS SNPs are in eQTL dataset (too rare)"

	return combined

if __name__ == "__main__":
	main()
