import sys
import subprocess
import pandas as pd
from shutil import copyfile
import os

(in_file,
	out_dir,
	out_file,
	tmp_raw,
	trait1_vcf,
	trait2_vcf,
	window,
	trait1_N,
	trait2_N,
	max_causal,
	verbose,
	save_finemap_threshold) = sys.argv[1:]

save_finemap_threshold = float(save_finemap_threshold)

trait1 = in_file.split("TRAIT1-")[1].split(".TRAIT2")[0]
trait2 = in_file.split("TRAIT2-")[1].split("/locus.")[0]

locus = in_file.split("locus.")[1].split(".txt")[0]
tmp_dir = f"{tmp_raw}/{locus}/"

def main():

	os.makedirs(out_dir, exist_ok=True)
	os.makedirs(tmp_dir, exist_ok=True)

	# Load data
	data = pd.read_csv(in_file, sep="\t")

	N_snps = data.shape[0]

	# Prep files for running finemap
	data = prep_finemap(data)

	# Check if there was an error in prep, like no remaining SNPS -- right now we don't really
	# do anything to log or handle it though
	if isinstance(data, str):
		return

	launch_finemap()

	min_pvalue_trait1 = min(list(data["pvalue_trait1"]))
	min_pvalue_trait2 = min(list(data["pvalue_trait2"]))
	feature1 = data["feature_trait1"].iloc[0]
	feature2 = data["feature_trait2"].iloc[0]
	chrom = data["chr_trait1"].iloc[0]
	snp_pos = data["seed_pos"].iloc[0]
	process_finemap_results(min_pvalue_trait1, min_pvalue_trait2, feature1, feature2, chrom, snp_pos, N_snps)

# Separates the preparation of data for finemap from the actual 
# launching of finemap. This is done to avoid duplicating code,
# since eCAVIAR and FINEMAP use very similar setups.
def prep_finemap(data):
	

	# Make required directories for FINEMAP eCAVIAR analysis

	os.makedirs(f"{tmp_dir}/vcftools", exist_ok=True)
	os.makedirs(f"{tmp_dir}/plink", exist_ok=True)
	os.makedirs(f"{tmp_dir}/ecaviar", exist_ok=True)
	os.makedirs(f"{tmp_dir}/finemap", exist_ok=True)

	# Get VCF file paths

	combined = data.copy()

	# Do we have one ref genome or two?
	two_refs = len([v for v in list(combined.columns.values) if "_vcf2" in v]) > 0

	# Two different cases depending on whether trait1 and trait2
	# are using same reference genome.
	if not two_refs:

		vcf1 = combined[[v for v in list(combined.columns.values) if "_vcf1" in v]]
		
		removal_list = compute_ld(vcf1, "_vcf1")

		if isinstance(removal_list, str):
			return removal_list

		# LD will just be the same for both
		copyfile(f"{tmp_dir}/ecaviar/data_vcf1.fixed.ld", f"{tmp_dir}/ecaviar/data_vcf2.fixed.ld")

		# Remove indices that produced NaNs in the LD computations
		removal_list = list(set(removal_list))
		combined = combined.drop(combined.index[removal_list])

	else:

		vcf1 = combined[[v for v in list(combined.colmns.values) if "_vcf1" in v]]
		vcf2 = combined[[v for v in list(combined.colmns.values) if "_vcf2" in v]]

		# Run PLINK on both VCFs.
		while True:
			removal_list = compute_ld(combined, "_vcf1")
			if isinstance(removal_list, str):
				return removal_list
			
			extension_list = compute_ld(combined, "_vcf2")
			if isinstance(extension_list, str):
				return extension_list
 
			removal_list.extend(extension_list)

			# Continue until no more NaNs
			if len(removal_list) == 0:
				break

			removal_list = list(set(removal_list))

			combined = combined.drop(combined.index[removal_list])

	# NOTE: it is possible that the above filters will make it so that there's
	# no longer a significant variant in one or both of the traits, but it will
	# run anyway unless a check for this is added

	# Reverse Z-scores for SNPs that are opposite of the VCF ref / alt designations
	
	flip_vcf1_indices = \
		(combined['effect_allele_trait1'] == combined['REF_vcf1']) & (combined['non_effect_allele_trait1'] == combined['ALT_vcf1'])
	if "REF_vcf2" in list(combined.columns.values):
		flip_vcf2_indices = \
			(combined['effect_allele_trait1'] == combined['REF_vcf2']) & (combined['non_effect_allele_trait1'] == combined['ALT_vcf2'])
	else:
		flip_vcf2_indices = flip_vcf1_indices

	
	combined = combined.reset_index(drop=True)
	# Apply it...

	if "REF_vcf2" in list(combined.columns.values):
		ref_vcf2 = "REF_vcf2"
	else:
		ref_vcf2 = "REF_vcf1"

	for i in range(combined.shape[0]):
	
		# We don't need to fix the whole matrix; we just need to make it so
		# the Z-scores being written for the trait are consistent in sign
		# with the LD matrix

		# We assume here that trait1 and trait2 alleles already align
		# because this was done in the harmonization step of preprocessing.

		# Breaking alignment between the two traits does not matter here 
		# since they're being fine-mapped independently here.

		# NOTE: Don't really need to realign the effect / non_effect
		# column for FINEMAP v1.4 -- see note below

		if combined['effect_allele_trait1'][i] == combined['REF_vcf1'][i]:
			combined.at[i, 'beta_trait1'] *= -1
	
		if combined['effect_allele_trait1'][i] == combined[ref_vcf2][i]:
			combined.at[i, 'beta_trait2'] *= -1
		
	# From the FINEMAP documentation for v1.4:
	# Columns beta and se are required for fine-mapping. Column maf is needed to output posterior 
	# effect size estimates on the allelic scale. All other columns are not required for computations and can be specified arbitrarily.

	# This is why we don't care that the SNP REF/ALT has not been sufficiently reversed above

	combined["ref_af_vcf1"] = combined["ref_af_vcf1"].apply(lambda x: min(x, 1-x))

	if "ref_af_vcf2" in list(combined.columns.values):
		ref_af_vcf2 = "ref_af_vcf2"
		combined["ref_af_vcf2"] = combined["ref_af_vcf2"].apply(lambda x: min(x, 1-x))
	else:
		ref_af_vcf2 = "ref_af_vcf1"

	# Then write Z-scores to file for FINEMAP
	with open(f"{tmp_dir}/ecaviar/trait1.z", "w") as w:
		snps = combined[['rsid', 'chr_trait1', 'snp_pos', 'effect_allele_trait1', 'non_effect_allele_trait1', 'ref_af_vcf1', 'beta_trait1', 'se_trait1']]
		snps.to_csv(w, index=False, header=['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se'], sep=" ")
	with open(f"{tmp_dir}/ecaviar/trait2.z", "w") as w:
		snps = combined[['rsid', 'chr_trait2', 'snp_pos', 'effect_allele_trait1', 'non_effect_allele_trait1', ref_af_vcf2, 'beta_trait2', 'se_trait2']]
		snps.to_csv(w, index=False, header=['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se'], sep=" ")
	return(combined)

# This function contains the code that's specific to FINEMAP,
# not shared with eCAVIAR.
def launch_finemap():

	# Write config file for finemap
	with open(f'{tmp_dir}/ecaviar/finemap.in', "w") as w:
		subprocess.run(f'echo z;ld;snp;config;cred;log;n_samples'.split(), stdout=w)

	with open(f'{tmp_dir}/ecaviar/finemap.in', "a") as a:
		subprocess.run(f'echo {tmp_dir}/ecaviar/trait1.z;{tmp_dir}/ecaviar/data_vcf1.fixed.ld;{tmp_dir}/ecaviar/finemap.trait1.snp;{tmp_dir}/ecaviar/finemap.trait1.config;{tmp_dir}/ecaviar/finemap.trait1.cred;{tmp_dir}/ecaviar/finemap.trait1.log;{trait1_N}'.split(), stdout=a)

	with open(f'{tmp_dir}/ecaviar/finemap.in', "a") as a:
		subprocess.run(f'echo {tmp_dir}/ecaviar/trait2.z;{tmp_dir}/ecaviar/data_vcf2.fixed.ld;{tmp_dir}/ecaviar/finemap.trait2.snp;{tmp_dir}/ecaviar/finemap.trait2.config;{tmp_dir}/ecaviar/finemap.trait2.cred;{tmp_dir}/ecaviar/finemap.trait2.log;{trait2_N}'.split(), stdout=a)

	# Run FINEMAP
	if verbose == True:
		subprocess.run(f'../bin/finemap/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --sss --in-files {tmp_dir}/ecaviar/finemap.in --n-causal-snps {max_causal} --n-iter 1000000 --n-conv-sss 1000'.split())
	else:
		subprocess.check_call(f'../bin/finemap/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --sss --in-files {tmp_dir}/ecaviar/finemap.in --n-causal-snps {max_causal} --n-iter 1000000 --n-conv-sss 1000'.split(), stdout=subprocess.DEVNULL)	
   
def process_finemap_results(min_pvalue_trait1, min_pvalue_trait2, feature1, feature2, chrom, snp_pos, N_snps):
	# Parse FINEMAP results to compute CLPP score
	trait1_probs = []
	trait2_probs = []
	with open(f"{tmp_dir}/ecaviar/finemap.trait1.snp") as f:
		f.readline()
		for line in f:
			data = line.strip().split()
			trait1_probs.append((int(data[0]), float(data[10])))
	with open(f"{tmp_dir}/ecaviar/finemap.trait2.snp") as f:
		f.readline()
		for line in f:
			data = line.strip().split()
			trait2_probs.append((int(data[0]), float(data[10])))
	
	trait1_probs = sorted(trait1_probs)
	trait2_probs = sorted(trait2_probs)

	# Make sure the SNPs have been paired up correctly
	assert len(trait1_probs) == len(trait2_probs)
	for i in range(len(trait1_probs)):
			assert trait1_probs[i][0] == trait2_probs[i][0]

	clpp = sum([trait1_probs[i][1]*trait2_probs[i][1] for i in range(len(trait1_probs))])
	# TODO: For multiple variants, it'll be something like this but it's not quite right...I think
	# the issue is that the individual variants are not mutually independent; one being causal affects the prob
	# of the others being causal.
	# clpp = 1 - reduce(mul, [1-(trait1_probs[i][1]*trait2_probs[i][1]) for i in range(len(trait1_probs))])

	# Note: LD for CLPP-mod score is being calculated only from the trait1 VCF;
	# may want to change this eventually to include the trait2 VCF too
	ld_file = f"{tmp_dir}/ecaviar/data_vcf1.fixed.ld"
	clpp_mod = get_clpp_mod(trait1_probs, trait2_probs, ld_file)

	# TODO: Figure the actual way to arrange the output files...
	if save_finemap_threshold != -1:
		
		if clpp_mod > save_finemap_threshold:
			os.makedirs(f"{out_dir}/finemap/finemap_probs/", exist_ok=True)

			copyfile(f"{tmp_dir}/ecaviar/finemap.trait1.snp", f"{out_dir}/finemap/finemap_probs/finemap.trait1.{locus}.snp")

			copyfile(f"{tmp_dir}/ecaviar/finemap.trait2.snp", f"{out_dir}/finemap/finemap_probs/finemap.trait2.{locus}.snp")

	# Write FINEMAP results to the desired file
	with open(out_file, "a") as a:
		a.write(f"{chrom}\t{snp_pos}\t{trait1}\t{trait2}\t{min_pvalue_trait1}\t{min_pvalue_trait2}\t{feature1}\t{feature2}\t{N_snps}\t{clpp}\t{clpp_mod}\n")

# Run PLINK on the locus of interest
def compute_ld(in_data, suffix):

	# We don't want to modify the input VCF within this function
	vcf = in_data.copy()

	# Repeatedly compute LD until we've eliminated all NaNs.
	removal_list = []
	while True:
		if vcf.shape[0] == 0:
			return "Fail: All SNPs have been eliminated through VCF filtering."

		sub_header = list(vcf.columns.values)

		# Write VCF to tmp file
		vcf.to_csv(f'{tmp_dir}/plink/data{suffix}.vcf', sep="\t", index=False, header=[h.replace(suffix, "") for h in sub_header])

		# Use PLINK to generate bim bam fam files
		command = f'''plink --vcf {tmp_dir}/plink/data{suffix}.vcf --keep-allele-order --make-bed --double-id --out {tmp_dir}/plink/data{suffix}_plinked'''.split()
		subprocess.run(command, stdout = subprocess.DEVNULL)

		# Use PLINK to generate LD score
		command = f'''plink -bfile {tmp_dir}/plink/data{suffix}_plinked --r square --out {tmp_dir}/ecaviar/data{suffix}'''.split()
		subprocess.run(command, stdout = subprocess.DEVNULL)

		# See if nans remain. If so, remove the offending lines.
		try:
			lines = [int(n.split(":")[0])-1 for n in subprocess.run(f"grep -n nan {tmp_dir}/ecaviar/data{suffix}.ld".split(), capture_output=True).stdout.decode('utf-8').strip().split("\n")]
		except:
			break

		# Save IDs of items being removed
		removal_list.extend(list(vcf.iloc[lines][f'rsid{trait1}']))
		
		# Remove desired rows (variants)
		vcf = vcf.drop(vcf.index[lines])

	# Get indices of items to remove in original list
	removal_list = [list(input_vcf[f'rsid{trait1}']).index(id) for id in removal_list]

	# Replace tabs with spaces because FINEMAP requires this.
	with open(f"{tmp_dir}/ecaviar/data{suffix}.fixed.ld", "w") as w:
		sp = subprocess.run(["cat", f"{tmp_dir}/ecaviar/data{suffix}.ld"], capture_output=True).stdout.decode("utf-8")
		w.write(sp.replace("\t", " "))

	return(removal_list)

# Compute LD-modified CLPP score.
# 
# The i_th variant of trait1_probs in this function must be the
# same as the i_th variant of trait2_probs.
def get_clpp_mod(trait1_probs, trait2_probs, ld_file):
	
	ld = []
	with open(ld_file) as f:
		for line in f:
			ld.append([float(f) for f in line.strip().split()]) 

	# Get modified CLPP score
	clpp_mod = 0
	for i in range(len(trait1_probs)):
		for j in range(len(trait2_probs)):
			snp_ld = ld[i][j]**2
			clpp_mod += trait1_probs[i][1] * trait2_probs[j][1] * snp_ld
	
	return clpp_mod


if __name__ == "__main__":
	main()
