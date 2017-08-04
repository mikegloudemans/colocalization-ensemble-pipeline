# Author: Mike Gloudemans
#
# Run FINEMAP algorithm for colocalization.
# Store full results in the output directory for this run.
# Return the CLPP score.
#

import subprocess

def run_finemap(locus, window=500000):
    
        conditional_level = 0   # Currently not used at all, but we may re-add this functionality later.

        gwas_chrom = locus.snp[0][3:]
        gwas_pos = locus.snp[1]
        combined = locus.data

	# Make required directories for eCAVIAR analysis
	subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix), shell=True)
	subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix), shell=True)
	subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix), shell=True)
	
	# For now, for simplicity, remove all positions that appear multiple times in the GWAS table.
	# This will avoid problems later in the pipeline, and doesn't remove too many SNPs anyway.
	dup_counts = {}
	for pos in combined['snp_pos']:
		dup_counts[pos] = dup_counts.get(pos, 0) + 1

	combined['dup_counts'] = [dup_counts[pos] for pos in combined['snp_pos']]
	combined = combined[combined['dup_counts'] == 1]

	snps = combined[['chr_y', 'snp_pos']]	
	# Write list of SNPs to a file for vcftools
	with open("/users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level) , "w") as w:
		snps.to_csv(w, index=False, header=False, sep="\t")

	# Get the region of interest from 1K genomes VCFs using tabix
	subprocess.check_call("tabix -h /mnt/lab_data/montgomery/shared/1KG/ALL.chr{1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz {1}:{6}-{7} > /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_prefiltered.recode_level{5}.vcf".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level, gwas_pos-window, gwas_pos+window), shell=True)

	# Use VCFtools to filter down to appropriate sites
	# (For the sake of a speedy analysis, we thought about requiring MAF above 0.01 in 1K Genomes, since hard to find eQTLs otherwise.
	# However, a visual check on this revealed that very few variants were removed at this filtering level,
	# possibly because these variants have also been filtered in GTEx.)
	command = 'vcftools --vcf /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_prefiltered.recode_level{5}.vcf --positions /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt --recode --recode-INFO-all --out /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes'.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level)
        subprocess.check_call(command, shell=True)

        # TODO: Move this to a separate function.
	# Loop through the output file, saving only the sites that appear in both the VCF
	# and in the combined SNPs list a single time (no more, no less!)
	used = set([])
	saved_list = set([])
	# Find SNPs that appear exactly once in the VCF
	with open('/users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.recode.vcf'.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level)) as f:
		for line in f:
			if line.startswith("#"):
				continue
			else:
				data = line.strip().split()
				pos = (int(data[0]), int(data[1]))
				saved_list.add(pos)
				if pos in used:
					saved_list.remove(pos)
				used.add(pos)

	# Remove SNPs from the VCF if they appear more than once
	with open('/users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.matched.recode.vcf'.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level), "w") as w:
		with open('/users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.recode.vcf'.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level)) as f:
			for line in f:
				if line.startswith("#"):
					w.write(line)
				else:
					data = line.strip().split()
					pos = (int(data[0]), int(data[1]))
					if pos in saved_list:
						w.write(line)

	# Remove rows from the SNP table if they don't appear in the VCF
	combined['pos_tuple'] = zip(combined['chr_y'], combined['snp_pos'])
	combined = combined[combined['pos_tuple'].isin(saved_list)]

	# Use PLINK to generate bim bam fam files
	command = '''/srv/persistent/bliu2/tools/plink_1.90_beta3_linux_x86_64/plink --vcf /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.matched.recode.vcf --keep-allele-order --make-bed --double-id --out /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes_plinked'''.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level)
	subprocess.check_call(command, shell=True)

	# Use PLINK to generate LD score
	command = '''/srv/persistent/bliu2/tools/plink_1.90_beta3_linux_x86_64/plink -bfile /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes_plinked --r square --out /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}'''.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level)
	subprocess.check_call(command, shell=True)

	# Fix LD-score by replacing nan values with 0.
	# TODO: Verify that this is valid and doesn't screw up results.
	# Also replace tabs with spaces because FINEMAP requires this.
	subprocess.check_call("sed s/nan/0/g /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.ld | sed s/\\\\t/\\ /g > /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level), shell=True)

	# Print z-scores to input files for FINEMAP.
	combined['ZSCORE_eqtl'] = combined['t-stat']
	# Figure out whether GWAS scores are in odds ratio or beta-se format
	if 'or' in combined:
		combined['ZSCORE_gwas'] = (combined['or']-1) / combined['se']
	else:
		combined['ZSCORE_gwas'] = (combined['beta_x']) / combined['se']
	with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level), "w") as w:
		snps = combined[['snp_pos', 'ZSCORE_eqtl']]
		snps.to_csv(w, index=False, header=False, sep=" ")

 	with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level), "w") as w:
		snps = combined[['snp_pos', 'ZSCORE_gwas']]	
		snps.to_csv(w, index=False, header=False, sep=" ")

	# Write config file for finemap
	# TODO TODO TODO TODO TODO: Right now we're just arbitrarily saying 5000 individuals in config just to get this running.
	# Need to do a bit more investigation into how the number of individuals used here affects
	# the results, and into how important it is to use the proper LD computations (e.g. GTEx-specific)
	# when exploring results. This is very important if we want to trust results.
	subprocess.check_call('echo "z;ld;snp;config;n-ind" > /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level), shell=True)
	subprocess.check_call('echo "/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.config;50000" >> /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level), shell=True)
	subprocess.check_call('echo "/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.config;500" >> /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level), shell=True)
	
	# Run FINEMAP
	subprocess.check_call('finemap --sss --in-files /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in --n-causal-max 1 --n-iterations 1000000 --n-convergence 50000'.format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level), shell=True)

	# Parse FINEMAP results to compute CLPP score
	gwas_probs = []
	eqtl_probs = []
	with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level)) as f:
		f.readline()
		for line in f:
			data = line.strip().split()
			gwas_probs.append((int(data[1]), float(data[2])))
        with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level)) as f:
		f.readline()
		for line in f:
			data = line.strip().split()
			eqtl_probs.append((int(data[1]), float(data[2])))
	gwas_probs = sorted(gwas_probs)
	eqtl_probs = sorted(eqtl_probs)

	assert len(gwas_probs) == len(eqtl_probs)
	for i in range(len(gwas_probs)):
		assert gwas_probs[i][0] == eqtl_probs[i][0]

	finemap_clpp = sum([gwas_probs[i][1] * eqtl_probs[i][1] for i in range(len(gwas_probs))])

	# Write FINEMAP results to the desired file
        # Note: appending will always work, since results always go to a different directory for each run.
        with open("{0}/{1}_finemap_clpp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
                a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene, conditional_level, snps.shape[0], finemap_clpp))

	# Remove temporary intermediate files to save space
	purge_tmp_files(locus, gwas_chrom, gwas_pos)

        return finemap_clpp

# Function purge_tmp_files
# Remove temporary files created during this run of eCAVIAR,
# to free up space.
def purge_tmp_files(locus, gwas_chrom, gwas_pos):
	
	# Why not just purge everything from the entire directory of the GWAS that we're working on?
	# Answer: it's possible that other jobs are still using those files.
	# Only purge the specific files created in this run.

        # (To avoid conflicts here, make sure each gwas/eqtl/snp triplet is run in its own process;
        # doubling up would cause problems).

	subprocess.call("rm -r /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3} 2> /dev/null".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene), shell=True)
	subprocess.call("rm -r /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3} 2> /dev/null".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene), shell=True)
	subprocess.call("rm -r /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3} 2> /dev/null".format(locus.gwas_suffix, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene), shell=True)
