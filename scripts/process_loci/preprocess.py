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
maf_ref_file = sys.argv[9]

def main():

	source_data = get_sumstats(source_file, chrom, pos, "_source", source_trait)
	lookup_data = get_sumstats(lookup_file, chrom, pos, "_lookup", lookup_trait)

	if isinstance(source_data, str) or source_data.shape[0] == 0:
		print(f"No overlapping data for {source_file} {lookup_file} {chrom} {pos} {source_trait} {lookup_trait}")
		return
	if isinstance(lookup_data, str) or lookup_data.shape[0] == 0:
		print(f"No overlapping data for {source_file} {lookup_file} {chrom} {pos} {source_trait} {lookup_trait}")
		return

	merge_data = combine_sumstats(source_data, lookup_data)
	if isinstance(merge_data, str) or merge_data.shape[0] == 0:
		print(f"No overlapping data for {source_file} {lookup_file} {chrom} {pos} {source_trait} {lookup_trait}")
		return

	if "effect_allele_source" in list(merge_data.columns.values) and "effect_allele_lookup" in list(merge_data.columns.values) and \
		"non_effect_allele_source" in list(merge_data.columns.values) and "non_effect_allele_lookup" in list(merge_data.columns.values):
		merge_data = harmonize_alleles(merge_data)

	merge_data = get_mafs(merge_data, maf_ref_file)

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
	combined = combined.drop(columns=['dup_counts'])

	# Check to make sure there are SNPs remaining; if not, just move on
	# to next gene.
	if combined.shape[0] == 0:
		return "No overlapping SNPs in eQTL and GWAS"

	# Check to make sure we still have significant GWAS hits and eQTLs, if desired
	#if "screening_thresholds" in settings and "gwas" in settings["screening_thresholds"]:
	#	if min(combined['pvalue_gwas']) > settings["screening_thresholds"]["gwas"]:
	#		return "No significant GWAS SNPs are in eQTL dataset (too rare)"

	return combined

# Re-align the allele directions of the source and lookup files so that the effect
# and non-effect alleles are consistent. Also, delete any entries for which the 
# alleles simply don't match.
def harmonize_alleles(merge_data):

	data = merge_data.copy().reset_index(drop=True)

	match_index = pd.Series([True]*data.shape[0])
	for i in data.index:
		if (data['effect_allele_source'][i] == data['effect_allele_lookup'][i] and
			data['non_effect_allele_source'][i] == data['non_effect_allele_lookup'][i]):
			# Everything agrees already
			pass

		elif (data['effect_allele_source'][i] == data['non_effect_allele_lookup'][i] and
						data['non_effect_allele_source'][i] == data['effect_allele_lookup'][i]):
			# Alleles are matching, but backwards

			tmp = data.at[i, 'effect_allele_source']
			data.at[i, 'effect_allele_source'] = data.at[i, 'non_effect_allele_source']
			data.at[i, 'non_effect_allele_source'] = tmp
	
			# Change OR, beta, and direction columns if needed
			if "beta_source" in list(data.columns.values):
				data.at[i, 'beta_source'] *= -1
			if "ZSCORE_source" in list(data.columns.values):
				data.at[i, 'ZSCORE_source'] *= -1
			if "tstat_source" in list(data.columns.values):
				data.at[i, 'tstat_source'] *= -1
			if "or_source" in list(data.columns.values):
				data.at[i, 'or_source'] = 1 / data.at[i, 'or_source']
			if "direction" in list(data.columns.values):
				def flip(x):
					if x == "-":
						return "+"
					if x == "+":
						return "-"
				data.at[i, 'direction'] = flip(data.at[i, 'direction'])

		else:
			# Not actually the same SNP; we will discard this position
			match_index.update(pd.Series([False], index=[i]))
	data = data[match_index]	
	data = data.drop(columns=['non_effect_allele_lookup', 'effect_allele_lookup'])	# Now they are the same as source anyway
	return data

# Load MAFs from a reference VCF, for every SNP in "data". For now, the direction
# is not important ; i.e. MAF = n and MAF = (1-n) are treated equivalently by
# COLOC.
def get_mafs(dataframe, maf_ref_file):

	chr_prefix = has_chr_prefix(maf_ref_file)
	if chr_prefix:
		chr_text = f"chr{chrom}"
	else:
		chr_text = chrom

	# Standard header is as follows, but let's make sure
	header = ["chrom_vcf", "pos_vcf", "rsid_vcf", "ref_vcf", "alt_vcf", "info_vcf"]
	with gzip.open(maf_ref_file, 'rb') as f:
		for line in f:
			if line.decode('utf-8').startswith("#CHROM"):
				data = line.decode('utf-8').strip().split()
				chrom_index = 0
				pos_index = data.index('POS')
				rsid_index = data.index('ID')
				ref_index = data.index('REF')
				alt_index = data.index('ALT')
				info_index = data.index('INFO')
				break

	stream = StringIO(subprocess.run(f"tabix {maf_ref_file} {chr_text}:{pos-window}-{pos+window} | cut -f{chrom_index+1},{pos_index+1},{rsid_index+1},{ref_index+1},{alt_index+1},{info_index+1}", capture_output=True, shell=True).stdout.decode('utf-8'))
		
	# For readability, load the header too
	# Load with pandas
	vcf = pd.read_csv(stream, sep="\t", names=header)

	# Remove variants not in the GWAS table
	vcf["pos_vcf"] = vcf["pos_vcf"].astype(int)

	vcf = vcf[vcf["pos_vcf"].isin(list(dataframe["snp_pos"]))]

	# Remove variants with position appearing multiple times
	dup_counts = {}
	for v in vcf["pos_vcf"]:
		dup_counts[v] = dup_counts.get(pos, 0) + 1
	vcf["dup_counts"] = [dup_counts[v] for v in vcf['pos_vcf']]
	vcf = vcf[vcf["dup_counts"] == 1]
	vcf = vcf.drop(columns=['dup_counts'])

	# Remove multiallelic variants with only one entry in VCF
	l = lambda x: "," not in x
	vcf = vcf[vcf["ref_vcf"].apply(l) & vcf["alt_vcf"].apply(l)]

	# Remove monoallelic variants.
	
	# Allele frequency might be input as counts or as percentages,
	# so handle this.
	min_af = 0.01
	example_info = list(vcf["info_vcf"])[0]
	if "AF" in example_info:
		def fn(x):
			info = [s for s in x.split(";") if s.startswith("AF=")][0]
			af = float(info.split("=")[1])
			return af
		vcf['ref_af'] = vcf["info_vcf"].apply(fn)
		vcf = vcf[(vcf['ref_af'] > min_af) & (1-vcf['ref_af'] > min_af)]
	elif "AC" in example_info and "N" in example_info:
		# Get allele count for minor allele, divided total number of alleles measured in population
		def fn(x):
			n_info = [s for s in x.split(";") if s.startswith("N=")][0]
			n = int(info.split("=")[1])
	
			info = [s for s in x.split(";") if s.startswith("AC=")][0]
			ac = int(info.split("=")[1])

			af = ac*1.0/n
			return af
		vcf['ref_af'] = vcf["info_vcf"].apply(fn)
		vcf = vcf[(vcf['ref_af'] > 0.01) & (1-vcf['ref_af'] > 0.01)]

	# Merge ref genome and MAFs with the sumstats dataframe frame
	merged = pd.merge(dataframe, vcf, left_on="snp_pos", right_on="pos_vcf")

	# Remove variants where alt/ref don't match between GWAS/eQTL and VCF
	# Flipped is okay. A/C and C/A are fine, A/C and A/G not fine.

	keep_indices = \
		(((merged['non_effect_allele_source'] == merged['ref_vcf']) & (merged['effect_allele_source'] == merged['alt_vcf'])) | \
		((merged['effect_allele_source'] == merged['ref_vcf']) & (merged['non_effect_allele_source'] == merged['alt_vcf']))) & \
		(((merged['non_effect_allele_lookup'] == merged['ref_vcf']) & (merged['effect_allele_lookup'] == merged['alt_vcf'])) | \
		((merged['effect_allele_lookup'] == merged['ref_vcf']) & (merged['non_effect_allele_lookup'] == merged['alt_vcf'])))

	merged = merged.drop(columns=['pos_vcf', 'rsid_vcf', 'alt_vcf', 'info_vcf'])

	merged = merged.reset_index(drop=True)
	merged = merged[keep_indices]

	return merged

# Helper function:
# Returns true if "chr" prefix in VCF
def has_chr_prefix(maf_ref_file):
	if maf_ref_file.endswith(".gz"):
		with gzip.open(maf_ref_file, 'rb') as vcf:
			for line in vcf:
				if line.decode('utf-8').startswith("#"):
					continue
				if line.decode('utf-8').startswith("chr"):
					chr_prefix = True
					break
				else:
					chr_prefix = False
					break
	else:
		with open(maf_ref_file, 'rb') as vcf:
			for line in vcf:
				if line.startswith("#"):
					continue
				if line.startswith("chr"):
					chr_prefix = True
					break
				else:
					chr_prefix = False
					break



if __name__ == "__main__":
	main()
