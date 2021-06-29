import sys
import subprocess
import pandas as pd
from shutil import copyfile

in_file, \
	out_dir,
	out_file,
	tmp_raw,
	source_vcf,
	lookup_vcf,
	window,
	source_N,
	lookup_N,
	max_causal,
	verbose,
	save_finemap_threshold = sys.argv[1:]

trait1 = in_file.split("TRAIT1-")[1].split(".TRAIT2")[0]
trait2 = in_file.split("TRAIT2-")[1].split("/locus.")[0]

locus = in_file.split("locus.")[1].split(".txt")
tmp_dir = f"{tmp_raw}/{locus}/"

def main():

	os.makedirs(out_dir)
	os.makedirs(tmp_dir)

	# Load data
	data = pd.read_csv(in_file, sep="\t", header=True)

	# Prep files for running finemap
	data = prep_finemap(data)

	launch_finemap(data)

	min_pvalue_source = min(list(data["source_pvalue"]))
	min_pvalue_lookup = min(list(data["lookup_pvalue"]))
	feature = data["feature"].iloc[0]

	process_finemap_results(min_pvalue_source, min_pvalue_lookup, feature)

# Separates the preparation of data for finemap from the actual 
# launching of finemap. This is done to avoid duplicating code,
# since eCAVIAR and FINEMAP use very similar setups.
def prep_finemap(data):
	

	# Make required directories for FINEMAP eCAVIAR analysis

	os.mkdirs(f"{tmp_dir}/vcftools")
	os.mkdirs(f"{tmp_dir}/plink")
	os.mkdirs(f"{tmp_dir}/ecaviar")
	os.mkdirs(f"{tmp_dir}/finemap")

	# Get VCF file paths

	combined = data.copy()

	# Do we have one ref genome or two?
	two_refs = len([v for v in list(combined.columns.values) if "_vcf2" in v]) > 0

	# Two different cases depending on whether source and lookup
	# are using same reference genome.
	if not two_refs:

		vcf1 = combined[[v for v in list(combined.columns.values) if "_vcf1" in v]]
		
		removal_list = compute_ld(vcf1, "_vcf1")

		if isinstance(removal_list, basestring):
			return removal_list

		# LD will just be the same for both
		copy_file(f"{tmp_dir}/ecaviar/lookup.fixed.ld", f"{tmp_dir}/ecaviar/source.fixed.ld")

		# Remove indices that produced NaNs in the LD computations
		removal_list = list(set(removal_list))
		combined = combined.drop(combined.index[removal_list])

	else:

		vcf1 = combined[[v for v in list(combined.colmns.values) if "_vcf1" in v]]
		vcf2 = combined[[v for v in list(combined.colmns.values) if "_vcf2" in v]]

		# Run PLINK on both VCFs.
		while True:
			removal_list = compute_ld(combined, "_vcf1")
			if isinstance(removal_list, basestring):
				return removal_list
			
			extension_list = compute_ld(combined, "_vcf2")
			if isinstance(extension_list, basestring):
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
		(combined['effect_allele_source'] == combined['ref_vcf1']) & (merged['non_effect_allele_source'] == merged['alt_vcf1'])
	if "ref_vcf2" in list(combined.columns.values):
		flip_vcf2_indices = \
			(combined['effect_allele_source'] == combined['ref_vcf2']) & (merged['non_effect_allele_source'] == merged['alt_vcf2'])
	else:
		flip_vcf2_indices = flip_vcf1_indices

	
	combined = combined.reset_index(drop=True)
	# Apply it...

	if "ref_vcf2" in list(combined.columns.values):
		ref_vcf2 = "ref_vcf2"
	else:
		ref_vcf2 = "ref_vcf1"

	for i in combined.shape[0]:
	
		# We don't need to fix the whole matrix; we just need to make it so
		# the Z-scores being written for the trait are consistent in sign
		# with the LD matrix

		if combined['effect_allele_source'][i] == combined['ref_vcf1'][i]:
			combined.at[i, 'ZSCORE_source'] *= -1
	
		if combined['effect_allele_lookup'][i] == combined[ref_vcf2][i]:
			combined.at[i, 'ZSCORE_lookup'] *= -1
		

	# Then write Z-scores to file for FINEMAP
	with open(f"{tmp_dir}/lookup.z", "w") as w:
		snps = combined[['snp_pos', 'ZSCORE_lookup']]
		snps.to_csv(w, index=False, header=False, sep=" ")
	with open("{tmp_dir}/source.z", "w") as w:
		snps = combined[['snp_pos', 'ZSCORE_source']]
		snps.to_csv(w, index=False, header=False, sep=" ")
	return(combined)

# This function contains the code that's specific to FINEMAP,
# not shared with eCAVIAR.
def launch_finemap():

	# Write config file for finemap
	subprocess.check_call(f'echo "z;ld;snp;config;n-ind" > {tmp_dir}/ecaviar/finemap.in')

	subprocess.check_call(f'echo "{tmp_dir}/ecaviar/source.z;
					{tmp_dir}/ecaviar/source.fixed.ld;{tmp_dir}/ecaviar/finemap.source.snp;{tmp_dir}/ecaviar/finemap.source.config;{source_N}" >> {tmp_dir}/ecaviar/finemap.in')

	subprocess.check_call(f'echo "{tmp_dir}/ecaviar/lookup.z;
					{tmp_dir}/ecaviar/lookup.fixed.ld;{tmp_dir}/ecaviar/finemap.lookup.snp;{tmp_dir}/ecaviar/finemap.lookup.config;{lookup_N}" >> {tmp_dir}/ecaviar/finemap.in')

	# Run FINEMAP
	if verbose == True:
		subprocess.check_call(f'../bin/finemap/finemap_v1.4_x86_64 --sss --in-files {tmp_dir}/ecaviar/finemap.in --n-causal-max {max_causal} --n-iterations 1000000 --n-convergence 1000')
	else:
		subprocess.check_call(f'../bin/finemap/finemap_v1.4_x86_64 --sss --in-files {tmp_dir}/ecaviar/finemap.in --n-causal-max {max_causal} --n-iterations 1000000 --n-convergence 1000 > /dev/null')	
   
def process_finemap_results(min_pvalue_source, min_pvalue_lookup, feature):
	# Parse FINEMAP results to compute CLPP score
	source_probs = []
	lookup_probs = []
	with open(f"{tmp_dir}/ecaviar/finemap.source.snp" as f:
		f.readline()
		for line in f:
			data = line.strip().split()
			source_probs.append((int(data[1]), float(data[2])))
	with open(f"{tmp_dir}/ecaviar/finemap.lookup.snp" as f:
		f.readline()
		for line in f:
			data = line.strip().split()
			lookup_probs.append((int(data[1]), float(data[2])))
	
	source_probs = sorted(source_probs)
	lookup_probs = sorted(lookup_probs)

	assert len(source_probs) == len(lookup_probs)
	for i in range(len(source_probs)):
			assert source_probs[i][0] == lookup_probs[i][0]

	clpp = 1 - reduce(mul, [1-(source_probs[i][1]*lookup_probs[i][1]) for i in range(len(source_probs))])
	ld_file = f"{tmp_dir}/ecaviar/source.fixed.ld"
	clpp_mod = get_clpp_mod(source_probs, lookup_probs, ld_file)

	# TODO: Figure the actual way to arrange the output files...
	if save_finemap_threshold != -1:
		
		if finemap_clpp_mod > save_finemap_threshold:
			os.makedirs(f"{out_dir}/finemap/finemap_probs/")

			copyfile(f"{tmp_dir}/ecaviar/finemap.source.snp", f"{out_dir}/finemap/finemap_probs/finemap.source.{locus}.snp")

			copyfile(f"{tmp_dir}/ecaviar/finemap.lookup.snp", f"{out_dir}/finemap/finemap_probs/finemap.lookup.{locus}.snp")

	# Write FINEMAP results to the desired file
	with open(out_file, "a") as a:
		a.write("{chrom}\t{snp_pos}\t{trait1}\t{trait2}\t{min_pvalue_source}\t{min_pvalue_lookup}\t{feature}\t{clpp}\t{clpp_mod}\n")

# Run PLINK on the locus of interest
def compute_ld(input_vcf, locus, suffix):

	# We don't want to modify the input VCF within this function
	vcf = input_vcf.copy()

	# Repeatedly compute LD until we've eliminated all NaNs.
	removal_list = []
	while True:
		if vcf.shape[0] == 0:
			return "Fail: All SNPs have been eliminated through VCF filtering."

		# Write VCF to tmp file
		vcf.to_csv(f'{tmp_dir}/plink/data{suffix}.vcf', sep="\t", index=False, header=True)

		# Use PLINK to generate bim bam fam files
		command = f'''plink2 --vcf {tmp_dir}/plink/data{suffix}.vcf --keep-allele-order --make-bed --double-id --out {tmp_dir}/plink/{0}/data{suffix}_plinked > /dev/null'''
		subprocess.check_call(command, shell=True)

		# Use PLINK to generate LD score
		command = f'''plink2 -bfile {tmp_dir}/plink/{0}/data{suffix}_plinked --r square --out {tmp_dir}/ecaviar/data{suffix} > /dev/null'''
		subprocess.check_call(command, shell=True)

		# See if nans remain. If so, remove the offending lines.
		try:
			lines = [int(n.split(":")[0])-1 for n in subprocess.check_output(f"grep -n nan {tmp_dir}/ecaviar/data{suffix}.ld", shell=True).strip().split("\n")]
		except:
			break

		# Save IDs of items being removed
		removal_list.extend(list(vcf.iloc[lines]['rsid{source}']))
		
		# Remove desired rows (variants)
		vcf = vcf.drop(vcf.index[lines])

	# Get indices of items to remove in original list
	removal_list = [list(input_vcf['rsid{source}']).index(id) for id in removal_list]

	# Replace tabs with spaces because FINEMAP requires this.
	subprocess.check_call(f"cat {tmp_dir}/ecaviar/data{suffix}.ld | sed s/\\\\t/\\ /g > {tmp_dir}/ecaviar/data{suffix}.fixed.ld", shell=True)

	return(removal_list)

# In this function, source and lookup probs should be pre-sorted
def get_clpp_mod(source_probs, lookup_probs, ld_file):
	
	ld = []
	with open(ld_file) as f:
		for line in f:
			ld.append([float(f) for f in line.strip().split()]) 


	for i in range(len(source_probs)):
		assert source_probs[i][0] == lookup_probs[i][0]

	# Get modified CLPP score
	clpp_mod = 0
	for i in range(len(source_probs)):
		for j in range(len(lookup_probs)):
			snp_ld = ld[i][j]**2
			clpp_mod += source_probs[i][1] * lookup_probs[j][1] * snp_ld
	
	return clpp_mod


if __name__ == "__main__":
	main()
