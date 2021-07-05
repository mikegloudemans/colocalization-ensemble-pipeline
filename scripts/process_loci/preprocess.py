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
trait1_file = sys.argv[3]
trait2_file = sys.argv[4]
chrom = sys.argv[5]
pos = int(sys.argv[6])
feature1 = sys.argv[7]
feature2 = sys.argv[8]
vcf_ref_file1 = sys.argv[9]
vcf_ref_file2 = sys.argv[10]

def main():

	trait1_data = get_sumstats(trait1_file, chrom, pos, "_trait1", feature1)
	trait2_data = get_sumstats(trait2_file, chrom, pos, "_trait2", feature2)

	print(f"Processing {trait1_file} {feature1} {trait2_file} {feature2} {chrom} {pos}")

	if isinstance(trait1_data, str) or trait1_data.shape[0] == 0:
		print(f"No overlapping data for {trait1_file} {trait2_file} {chrom} {pos} {feature1} {feature2}")
		return
	if isinstance(trait2_data, str) or trait2_data.shape[0] == 0:
		print(f"No overlapping data for {trait1_file} {trait2_file} {chrom} {pos} {feature1} {feature2}")
		return

	merge_data = combine_sumstats(trait1_data, trait2_data)
	if isinstance(merge_data, str) or merge_data.shape[0] == 0:
		print(f"No overlapping data for {trait1_file} {trait2_file} {chrom} {pos} {feature1} {feature2}")
		return

	if "effect_allele_trait1" in list(merge_data.columns.values) and "effect_allele_trait2" in list(merge_data.columns.values) and \
		"non_effect_allele_trait1" in list(merge_data.columns.values) and "non_effect_allele_trait2" in list(merge_data.columns.values):
		merge_data = harmonize_alleles(merge_data)

	merge_data = get_ref_vcf(merge_data, vcf_ref_file1, "_vcf1")
	if vcf_ref_file1 != vcf_ref_file2:
		merge_data = get_ref_vcf(merge_data, vcf_ref_file2, "_vcf2")

	merge_data["seed_pos"] = pos 

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
			table = table.rename(index=str, columns={"trait": "feature"})
		elif "feature" in list(table.columns.values):
			table = table[table['feature'] == trait].copy()
		elif "gene" in list(table.columns.values):
			table = table[table['gene'] == trait].copy()
			table = table.rename(index=str, columns={"gene": "feature"})
		else:
			return("No trait column specified; trait filtering is impossible")
	else:
		table['feature'] = "none"
		
	if table.shape[0] == 0:
		return "No summary statistics found at this locus."

	table['pvalue'] = table['pvalue'].astype(float)

	if "ref" in list(table.columns.values):
		table = table.rename(index=str, columns={"ref": "non_effect_allele"})
		table = table.rename(index=str, columns={"alt": "effect_allele"})

	if "effect_direction" in list(table.columns.values):
		table = table.rename(index=str, columns={"effect_direction": "direction"})
				
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
		table = table[~pd.isna(table["se"])] # Thanks GTEx for making me have to do this
		table = table.rename(index=str, columns={"se": "se" + suffix})

	print(table.head(10))

	table['snp_pos'] = table['snp_pos'].astype(int)
	
	if min(table['pvalue']) < 1e-150:
		table['pvalue'] = table['pvalue'].apply(lambda x: max(x, 1e-150))

	# Sometimes the GWAS SNP is outside of the range of eQTLs tested for a certain
	# gene, or on the outside fringe of the range. If this is the case, then skip it.
	if pos > max(table['snp_pos']) - 50000 or pos < min(table['snp_pos']) + 50000:
		return "Sumstats interval too small."
	
	# We'll want to know the original "seed" position when logging results

	return table

def combine_sumstats(trait1_data, trait2_data):

	combined = pd.merge(trait1_data, trait2_data, on="snp_pos", suffixes=("_trait1", "_trait2"))
   
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

# Re-align the allele directions of the trait1 and trait2 files so that the effect
# and non-effect alleles are consistent. Also, delete any entries for which the 
# alleles simply don't match.
def harmonize_alleles(merge_data):

	data = merge_data.copy().reset_index(drop=True)

	match_index = pd.Series([True]*data.shape[0])
	for i in range(data.shape[0]):
		if (data['effect_allele_trait1'][i] == data['effect_allele_trait2'][i] and
			data['non_effect_allele_trait1'][i] == data['non_effect_allele_trait2'][i]):
			# Everything agrees already
			pass

		elif (data['effect_allele_trait1'][i] == data['non_effect_allele_trait2'][i] and
						data['non_effect_allele_trait1'][i] == data['effect_allele_trait2'][i]):
			# Alleles are matching, but backwards

			tmp = data.at[i, 'effect_allele_trait1']
			data.at[i, 'effect_allele_trait1'] = data.at[i, 'non_effect_allele_trait1']
			data.at[i, 'non_effect_allele_trait1'] = tmp
	
			# Change OR, beta, and direction columns if needed
			if "beta_trait1" in list(data.columns.values):
				data.at[i, 'beta_trait1'] *= -1
			if "ZSCORE_trait1" in list(data.columns.values):
				data.at[i, 'ZSCORE_trait1'] *= -1
			if "tstat_trait1" in list(data.columns.values):
				data.at[i, 'tstat_trait1'] *= -1
			if "or_trait1" in list(data.columns.values):
				data.at[i, 'or_trait1'] = 1 / data.at[i, 'or_trait1']
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
	data = data.drop(columns=['non_effect_allele_trait2', 'effect_allele_trait2'])	# Now they are the same as trait1 anyway
	return data

# Load MAFs and alleles from a reference VCF, for every SNP in "data". For now, the direction
# is not important ; i.e. MAF = n and MAF = (1-n) are treated equivalently by
# COLOC.
def get_ref_vcf(dataframe, vcf_ref_file, suffix):

	chr_prefix = has_chr_prefix(vcf_ref_file)
	if chr_prefix:
		chr_text = f"chr{chrom}"
	else:
		chr_text = chrom

	# Get header

	with gzip.open(vcf_ref_file, 'rb') as f:
		for line in f:
			if line.decode('utf-8').startswith("#CHROM"):
				header = line.decode('utf-8').strip().split()

				header = [h + suffix for h in header]
				break


	stream = StringIO(subprocess.run(f"tabix {vcf_ref_file} {chr_text}:{pos-window}-{pos+window}".split(), capture_output=True).stdout.decode('utf-8'))	

	# For readability, load the header too
	# Load with pandas
	vcf = pd.read_csv(stream, sep="\t", names=header)

	# Remove variants not in the GWAS table
	vcf[f"POS{suffix}"] = vcf[f"POS{suffix}"].astype(int)

	vcf = vcf[vcf[f"POS{suffix}"].isin(list(dataframe["snp_pos"]))]

	# Remove variants with position appearing multiple times
	dup_counts = {}
	for v in vcf[f"POS{suffix}"]:
		dup_counts[v] = dup_counts.get(pos, 0) + 1
	vcf["dup_counts"] = [dup_counts[v] for v in vcf[f'POS{suffix}']]
	vcf = vcf[vcf["dup_counts"] == 1]
	vcf = vcf.drop(columns=['dup_counts'])

	# Remove multiallelic variants with only one entry in VCF
	l = lambda x: "," not in x
	vcf = vcf[vcf[f"REF{suffix}"].apply(l) & vcf[f"ALT{suffix}"].apply(l)]

	# Remove monoallelic variants.
	
	# Allele frequency might be input as counts or as percentages,
	# so handle this.
	min_af = 0.01
	example_info = list(vcf[f"INFO{suffix}"])[0]
	if "AF" in example_info:
		def fn(x):
			info = [s for s in x.split(";") if s.startswith("AF=")][0]
			af = float(info.split("=")[1])
			return af
		vcf[f'ref_af{suffix}'] = vcf[f"INFO{suffix}"].apply(fn)
		vcf = vcf[(vcf[f'ref_af{suffix}'] > min_af) & (1-vcf[f'ref_af{suffix}'] > min_af)]
	elif "AC" in example_info and "N" in example_info:
		# Get allele count for minor allele, divided total number of alleles measured in population
		def fn(x):
			n_info = [s for s in x.split(";") if s.startswith("N=")][0]
			n = int(info.split("=")[1])
	
			info = [s for s in x.split(";") if s.startswith("AC=")][0]
			ac = int(info.split("=")[1])

			af = ac*1.0/n
			return af
		vcf[f'ref_af{suffix}'] = vcf[f"INFO{suffix}"].apply(fn)

		vcf = vcf[(vcf[f'ref_af{suffix}'] > min_af) & (1-vcf[f'ref_af{suffix}'] > min_af)]
	
	# Merge ref genome and MAFs with the sumstats dataframe frame
	merged = pd.merge(dataframe, vcf, left_on="snp_pos", right_on=f"POS{suffix}")

	# Remove variants where alt/ref don't match between GWAS/eQTL and VCF
	# Flipped is okay. A/C and C/A are fine, A/C and A/G not fine.
	
	# Only need to pay attention to the "trait1" file since they will have already been harmonized by now for "trait2" file
	keep_indices = \
		(((merged['non_effect_allele_trait1'] == merged[f'REF{suffix}']) & (merged['effect_allele_trait1'] == merged[f'ALT{suffix}'])) | \
		((merged['effect_allele_trait1'] == merged[f'REF{suffix}']) & (merged['non_effect_allele_trait1'] == merged[f'ALT{suffix}']))) 

	# ^ Probably need to harmonize these too for FINEMAP, but we haven't gotten there quite yet

	# NOTE: Might be useful to write out a "short" version of this table for steps like COLOC that don't
	# need all the genotype data, along with a "long" version for those like FINEMAP that do need it
	merged = merged.reset_index(drop=True)
	merged = merged[keep_indices]

	return merged

# Helper function:
# Returns true if "chr" prefix in VCF
def has_chr_prefix(ref_file):
	with gzip.open(ref_file, 'rb') as vcf:
		for line in vcf:
			if line.decode('utf-8').startswith("#"):
				continue
			if line.decode('utf-8').startswith("chr"):
				chr_prefix = True
				break
			else:
				chr_prefix = False
				break

if __name__ == "__main__":
	main()
