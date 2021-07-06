# Author: Abhiram Rao
# Author: Mike Gloudemans

import sys
import subprocess
import pandas as pd
import os

(in_file,
	out_file,
	tmp_raw,
	trait1_vcf,
	window,
	verbose) = sys.argv[1:]

trait1 = in_file.split("TRAIT1-")[1].split(".TRAIT2")[0]
trait2 = in_file.split("TRAIT2-")[1].split("/locus.")[0]

locus = in_file.split("locus.")[1].split(".txt")[0]
tmp_dir = f"{tmp_raw}/{locus}/"

def main():

	data = pd.read_csv(in_file, sep="\t")

	run_smr(data)

	parse_smr_results(data)

# Integration of code written by Ram
def run_smr(in_data):
	
	os.makedirs(f"{tmp_dir}/smr", exist_ok=True)
	
	data = in_data.copy()

	data = harmonize_mafs(data)

	# NOTE: At this point we're assuming that rsids, betas, ses, and ref_afs are all
	# already present in the data. This should probably be verified at an earlier
	# point

	smr_exe = '../bin/smr/smr_Linux'

	# For LD files, I'm going to have to generate them on the fly

	trait2_pval_threshold = 1 # p-value threshold so that we run SMR at every test
	maf_threshold = 0.01
	cis_wind = 10000 # large cis window to include all SNPs
	diff_freq = 1 # allowable allele freq difference between datasets
	diff_freq_prop = 1 # all SNPs are allowed to have allele freq differences of at least diff.freq

	########################
	# Run SMR and get pvalue
	########################

	# inferred info
	trait2 = list(data['feature_trait2'])[0]
	chrom = list(data["chr_trait1"])[0]
	ref_snp = list(data["seed_pos"])[0]
	
	# Create ESD, which contains information on trait 2, or the potential mediating
	# trait within the SMR model
	esd = pd.DataFrame({'Chr': chrom, 'SNP': data["rsid"], 'Bp': data["snp_pos"], 'A1': data["non_effect_allele_trait1"], 'A2': data["non_effect_allele_trait1"], 'Freq': data['ref_af_vcf1'], 'Beta': data["beta_trait2"], 'se': data["se_trait2"], 'p': data["pvalue_trait2"]})
	esd = esd[['Chr', 'SNP', 'Bp', 'A1', 'A2', 'Freq', 'Beta', 'se', 'p']]
	esd.to_csv(f'{tmp_dir}/smr/smr.esd', sep="\t", index = False)

	## create flist file
	# from the documentation: GeneticDistance does not matter
	# orientation also doesn't matter since we're not displaying it visually
	flist = pd.DataFrame({'Chr': [chrom], 'ProbeID': [trait2], 'GeneticDistance': ['0'], 'ProbeBp': [ref_snp], 'Gene': [trait2], 'Orientation': ['+'], 'PathOfEsd': [f'{tmp_dir}/smr/smr.esd']})
	flist.to_csv(f'{tmp_dir}/smr/smr.flist', sep="\t", index = False)

	## run SMR to create binary files
	subprocess.run(f'{smr_exe} --eqtl-flist {tmp_dir}/smr/smr.flist --make-besd --out {tmp_dir}/smr/smr.besd'.split())

	## run SMR with corresponding gwas file 
	data['n'] = "NA"	# It says in the SMR manual that this value doesn't even matter

	# Again, check to make sure the sign is correct
	gwas = data[['rsid','effect_allele_trait1','non_effect_allele_trait1','ref_af_vcf1','beta_trait1','se_trait1','pvalue_trait1','n']]
	gwas.columns = ['SNP','A1','A2','freq','b','se','p','n']
	gwas.to_csv(f"{tmp_dir}/smr/smr.gwas", sep="\t", index = False)

	# Run Plink to generate BIM BAM FAM files for the reference genome	
	make_bim_bam_fam(data)

	subprocess.run(f'{smr_exe} --bfile {tmp_dir}/smr/plink --gwas-summary {tmp_dir}/smr/smr.gwas --beqtl-summary {tmp_dir}/smr/smr.besd --maf {maf_threshold} --thread-num 1 --diff-freq {diff_freq} --diff-freq-prop {diff_freq_prop} --peqtl-smr {trait2_pval_threshold} --cis-wind {cis_wind} --out {tmp_dir}/smr/smr.out'.split())


def parse_smr_results(data):
	# Now parse results
	with open(f'{tmp_dir}/smr/smr.out.smr') as f:
		f.readline()
		data = f.readline().strip().split()
		smr_pval = float(data[18])
		smr_heidi_pval = float(data[19])

	# Get some metadata for the output file
	min_pvalue_trait1 = min(list(data["pvalue_trait1"]))
	min_pvalue_trait2 = min(list(data["pvalue_trait2"]))
	feature1 = data["feature_trait1"].iloc[0]
	feature2 = data["feature_trait2"].iloc[0]
	chrom = data["chr_trait1"].iloc[0]
	snp_pos = data["seed_pos"].iloc[0]

	N_snps = data.shape[0]

	# Write results to the output file
	with open(out_file, "a") as a:
		a.write(f"{chrom}\t{ref_snp}\t{trait1}\t{trait2}\t{min_pvalue_trait1}\t{min_pvalue_trait2}\t{feature1}\t{feature2}\t{N_snps}\t{smr_pval}\t{smr_heidi_pval}\n")

def make_bim_bam_fam(in_data):

	# We don't want to modify the input VCF within this function
	vcf = in_data.copy()
	vcf = vcf[[v for v in list(vcf.columns.values) if "_vcf1" in v]]
	vcf = vcf.drop(columns = ["ref_af_vcf1"])

	sub_header = list(vcf.columns.values)

	# Write VCF to tmp file
	vcf.to_csv(f'{tmp_dir}/smr/preplink.vcf', sep="\t", index=False, header=[h.replace("_vcf1", "") for h in sub_header])

	# Use PLINK to generate bim bam fam files
	command = f'''plink --vcf {tmp_dir}/smr/preplink.vcf --keep-allele-order --make-bed --double-id --out {tmp_dir}/smr/plink'''.split()
	subprocess.run(command, stdout = subprocess.DEVNULL)

def harmonize_mafs(data):
	
	# Make sure that the ref_af is referring to A1 and not A2 in the trait association sumstats,
	# since that's where it matters

	combined = data.copy()

	combined = combined.reset_index(drop=True)

	for i in range(combined.shape[0]):

		if combined['effect_allele_trait1'][i] == combined['REF_vcf1'][i]:
			combined.at[i, 'ref_af_vcf1'] = 1-combined.at[i, 'ref_af_vcf1']

	return(combined)

if __name__ == "__main__":
	main()
