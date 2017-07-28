#!/usr/bin/python
# Author: Mike Gloudemans
# Date: 10/10/2016

import os
import subprocess
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab
import math
import sys
import multiprocessing
import gzip
import collections

if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

# TODO TODO TODO: Think about how to make this much more object-oriented.
# How would I do this if starting from scratch? Shouldn't be too hard really.

# TODO: Add ArgumentParser to parse all of this input directly, in a separate function
# But actually TODO...just make this a module so we don't have to pass the arguments back
# and forth via command line at all.
window = 500000

# Read command line arguments
gwas_chrom = int(sys.argv[1].replace('chr', ''))
gwas_pos = int(sys.argv[2])
gwas_file = sys.argv[3]

gwas_prefix = "/".join(gwas_file.split("/")[:-1])
gwas_suffix = gwas_file.split("/")[-1].replace(".", "_")

gwas_threshold = float(sys.argv[4])

# Determine output directory
threshold_directory = '{0:.1e}'.format(gwas_threshold)
base_output_dir = "/users/mgloud/projects/brain_gwas/output/{0}/{1}".format(gwas_suffix.replace(".", "_"), threshold_directory)


def main():

	# Subset GWAS list to SNPs near the GWAS position
	gwas_table = pd.read_csv(gwas_file, sep="\t")
	gwas_table = gwas_table[(gwas_table['snp_pos'] > gwas_pos - window) & (gwas_table['snp_pos'] < gwas_pos + window)]
	gwas_table = gwas_table[(gwas_table['chr'] == gwas_chrom) | (gwas_table['chr'] == 'chr{0}'.format(gwas_chrom))]

	# Create base directory for output files, in case not already created
	subprocess.call("mkdir -p {0}".format(base_output_dir), shell=True)

        # TODO: make this generalize to an arbitrary list of input eqtl files.
        # Just do this analysis for every one of the input files.
        # NOTE: All the stuff in this section can be done just once for each locus,
        # and then all methods run on the loaded data.
	for tissue in os.listdir("/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_data/MatrixEQTL/allCisSNPGenePairs/"):

		if "Brain" not in tissue:
			break

		tissue_prefix = "_".join(tissue.split(".")[0].split("_")[:-1])
		
		# Get eQTL data using tabix
		header = subprocess.check_output("zcat /users/mgloud/projects/brain_gwas/data/" + \
			"eqtls/tabix/{0}_Analysis.cis.eqtl.gz | head -n 1".format(tissue_prefix, \
			 gwas_chrom, gwas_pos - window, gwas_pos + window), shell=True)
		raw_eqtls = subprocess.check_output("tabix /users/mgloud/projects/brain_gwas/" + \
			 "data/eqtls/tabix/{0}_Analysis.cis.eqtl.gz -h {1}:{2}-{3}".format(tissue_prefix, \
			 gwas_chrom, gwas_pos - window, gwas_pos + window), shell=True)
		eqtls = pd.read_csv(StringIO(header + raw_eqtls), sep="\t")

		# Get full list of all genes whose eQTLs we're testing at this site
		genes = set(eqtls['gene'])

		# Run colocalization pipeline for every gene in the list.
                for gene in genes:
			gwas_eqtl_colocalization(gwas_chrom, gwas_pos, gene, eqtls, gwas_table, gwas_suffix, tissue_prefix)
		
### Function gwas_eqtl_colocalization
# Test whether a gwas SNP and an eQTL colocalize. Currently will
# always return True unless it's impossible to perform the test.
# In future, might return True/False to signal whether we should
# run fastQTL again conditioned on the top remaining SNP.
#
# Input: gwas location, gene to test for eQTL colocalization, set of eqtls
# Output: True if test ran, False if test did not occur.
def gwas_eqtl_colocalization(gwas_chrom, gwas_pos, gene, eqtls, gwas_table, gwas_suffix, \
	tissue_prefix):

	# Get the set of genes we're interested in
	eqtl_subset = eqtls[eqtls['gene'] == gene]
	
	if eqtl_subset.shape[0] == 0:
		print "EQTL_SUBSET_EMPTY_ERROR"
		print gwas_chrom, gwas_pos, gene, current_level, tissue_prefix
		return False

	# Sometimes the GWAS SNP is outside of the range of eQTLs tested for a certain
	# gene. If this is the case, then skip it.
	if gwas_pos > max(eqtl_subset['snp_pos']) + 10000 or gwas_pos < min(eqtl_subset['snp_pos'] - 10000):
		return False
	
	# Join the list of eQTL SNPs with the list of GWAS SNPs
	combined = pd.merge(gwas_table, eqtl_subset, on="snp_pos")

	# Check to make sure there are SNPs remaining; if not, just move on
	# to next gene.
	if combined.shape[0] == 0: 
		return False

	# Run eCAVIAR or FINEMAP on SNPs.
	run_finemap(combined, gwas_chrom, gwas_pos, tissue_prefix, gene)


# Function plot_colocalization_test
#
# Produces a paired Manhattan plot and a colocalization
# scatterplot for the given locus. 
#
# TODO: Go through this section line-by-line and see if it's still worth using.
# Something makes me think though, that this might best be a completely different function.
def plot_colocalization_test(gwas_chrom, gwas_pos, gene, combined_table, gwas_suffix, \
        tissue_prefix, current_level, clpp):

	subprocess.call("mkdir -p {0}/plots/{1}_{2}/{3}".format(base_output_dir, gwas_chrom, gwas_pos, tissue_prefix), shell=True)

	plt.figure(figsize=(10,10))
	plt.scatter([-1 * math.log10(p) for p in combined_table['pvalue_x']], [-1 * math.log10(p) for p in combined_table['pvalue_y']], c=combined_table['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)

	if max([-1 * math.log10(p) for p in combined_table['pvalue_x']] + [-1 * math.log10(p) for p in combined_table['pvalue_y']]) < 20:
		plt.axis([0, 20, 0, 20])
	else:
		plt.axis([0, max([20] + [-1 * math.log10(p) for p in combined_table['pvalue_x']] + [-1 * math.log10(p) for p in combined_table['pvalue_y']]), 0, max([20] + [-1 * math.log10(p) for p in combined_table['pvalue_x']] + [-1 * math.log10(p) for p in combined_table['pvalue_y']])])
	plt.xlabel('GWAS -log p-value', fontsize=16)
	plt.ylabel('eQTL -log p-value', fontsize=16)
	plt.title('{0} CLPP = {1}'.format(gene, clpp), fontsize=24)
	plt.savefig("{0}/plots/{1}_{2}/{3}/{4}.png".format(base_output_dir, gwas_chrom, gwas_pos, tissue_prefix, gene))
	plt.gcf().clear()
	plt.close()

	# Also create a LocusZoom-style plot showing the GWAS and eQTL signals next to one another.
	plt.figure(figsize=(20,10))
	plt.subplot(211)
	plt.scatter(combined_table['snp_pos'], [-1 * math.log10(p) for p in combined_table['pvalue_x']], c=combined_table['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)
	plt.plot((gwas_pos, gwas_pos), (-2, max([-1 * math.log10(p) for p in combined_table['pvalue_x']])), 'k--')
	plt.ylabel('GWAS -log p-value', fontsize=16)
	plt.subplot(212)
	plt.scatter(combined_table['snp_pos'], [-1 * math.log10(p) for p in combined_table['pvalue_y']], c=combined_table['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)
	plt.plot((gwas_pos, gwas_pos), (-2, max([-1 * math.log10(p) for p in combined_table['pvalue_y']])), 'k--')
	plt.ylabel('eQTL -log p-value', fontsize=16)
	plt.xlabel('Position', fontsize=16)
	plt.savefig("{0}/plots/{1}_{2}/{3}/{4}_manhattan.png".format(base_output_dir, gwas_chrom, gwas_pos, tissue_prefix, gene))
	plt.gcf().clear()
	plt.close()
	
def run_finemap(combined, gwas_chrom, gwas_pos, tissue_prefix, gene):

	# Make required directories for eCAVIAR analysis
	subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), shell=True)
	subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), shell=True)
	subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), shell=True)
	
        # Purge old eCAVIAR results
        # TODO: In the modified version, every results folder will be its own, so there
        # won't be any need to do this step.
	# This step is critical. If we don't call it, then eCAVIAR will append results to the 
	# previously generated list, which causes problems in later parsing.
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/*".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), shell=True)
  
	# For now, for simplicity, remove all positions that appear multiple times in the GWAS table.
	# This will avoid problems later in the pipeline, and doesn't remove too many SNPs anyway.
	dup_counts = {}
	for pos in combined['snp_pos']:
		dup_counts[pos] = dup_counts.get(pos, 0) + 1

	combined['dup_counts'] = [dup_counts[pos] for pos in combined['snp_pos']]
	combined = combined[combined['dup_counts'] == 1]

	snps = combined[['chr_y', 'snp_pos']]	
	# Write list of SNPs to a file for vcftools
	with open("/users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level) , "w") as w:
		snps.to_csv(w, index=False, header=False, sep="\t")

	# Get the region of interest from 1K genomes VCFs using tabix
	subprocess.check_call("tabix -h /mnt/lab_data/montgomery/shared/1KG/ALL.chr{1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz {1}:{6}-{7} > /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_prefiltered.recode_level{5}.vcf".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level, gwas_pos-window, gwas_pos+window), shell=True)

	# Use VCFtools to filter down to appropriate sites
	# (For the sake of a speedy analysis, we thought about requiring MAF above 0.01 in 1K Genomes, since hard to find eQTLs otherwise.
	# However, a visual check on this revealed that very few variants were removed at this filtering level,
	# possibly because these variants have also been filtered in GTEx.)
	command = 'vcftools --vcf /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_prefiltered.recode_level{5}.vcf --positions /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt --recode --recode-INFO-all --out /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes'.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)
        subprocess.check_call(command, shell=True)

	# TODO: Consider whether we should be using GTEx VCFs to compute LD instead of
	# using 1000 Genomes VCFs. Would have to compute LD separately for each tissue in this case

        # TODO: Move this to a separate function.
	# Loop through the output file, saving only the sites that appear in both the VCF
	# and in the combined SNPs list a single time (no more, no less!)
	used = set([])
	saved_list = set([])

	# Find SNPs that appear exactly once in the VCF
	with open('/users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.recode.vcf'.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)) as f:
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
	with open('/users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.matched.recode.vcf'.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), "w") as w:
		with open('/users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.recode.vcf'.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)) as f:
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
	command = '''/srv/persistent/bliu2/tools/plink_1.90_beta3_linux_x86_64/plink --vcf /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.matched.recode.vcf --keep-allele-order --make-bed --double-id --out /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes_plinked'''.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)
	subprocess.check_call(command, shell=True)

	# Use PLINK to generate LD score
	command = '''/srv/persistent/bliu2/tools/plink_1.90_beta3_linux_x86_64/plink -bfile /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes_plinked --r square --out /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}'''.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)
	subprocess.check_call(command, shell=True)

	# Fix LD-score by replacing nan values with 0.
	# TODO: Verify that this is valid and doesn't screw up results.
	# Also replace tabs with spaces because FINEMAP requires this.
	subprocess.check_call("sed s/nan/0/g /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.ld | sed s/\\\\t/\\ /g > /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)

	# Print z-scores to input files for FINEMAP.
	combined['ZSCORE_eqtl'] = combined['t-stat']
	# Figure out whether GWAS scores are in odds ratio or beta-se format
	if 'or' in combined:
		combined['ZSCORE_gwas'] = (combined['or']-1) / combined['se']
	else:
		combined['ZSCORE_gwas'] = (combined['beta_x']) / combined['se']
	with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), "w") as w:
		snps = combined[['snp_pos', 'ZSCORE_eqtl']]
		snps.to_csv(w, index=False, header=False, sep=" ")

 	with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), "w") as w:
		snps = combined[['snp_pos', 'ZSCORE_gwas']]	
		snps.to_csv(w, index=False, header=False, sep=" ")

	# NEW: FINEMAP pipeline

	# Write config file for finemap
	# TODO TODO TODO TODO TODO: Right now we're just arbitrarily saying 5000 individuals in config just to get this running.
	# Need to do a bit more investigation into how the number of individuals used here affects
	# the results, and into how important it is to use the proper LD computations (e.g. GTEx-specific)
	# when exploring results. This is very important if we want to trust results.
	subprocess.check_call('echo "z;ld;snp;config;n-ind" > /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.check_call('echo "/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.config;50000" >> /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.check_call('echo "/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.config;500" >> /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	
	# Run FINEMAP
	subprocess.check_call('finemap --sss --in-files /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in --n-causal-max 1 --n-iterations 1000000 --n-convergence 50000'.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)

	# Parse FINEMAP results to compute CLPP score
	gwas_probs = []
	eqtl_probs = []
	with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)) as f:
		f.readline()
		for line in f:
			data = line.strip().split()
			gwas_probs.append((int(data[1]), float(data[2])))
        with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)) as f:
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
        with open("{0}/{1}_finemap_clpp_status.txt".format(base_output_dir, gwas_suffix.replace(".", "_")), "a") as a:
                a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level, snps.shape[0], finemap_clpp))

	# Remove temporary intermediate files to save space
	purge_tmp_files(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)

	# Plot results
	plot_colocalization_test(gwas_chrom, gwas_pos, gene, combined, gwas_suffix, \
        	tissue_prefix, clpp)

### Function purge_tmp_files
# Remove temporary files created during this run of eCAVIAR,
# to free up space.
def purge_tmp_files(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level):
	
	# Why not just purge everything from the entire directory of the GWAS that we're working on?
	# Answer: it's possible that other jobs are still using those files.
	# Only purge the specific files created in this run.

	# Remove intermediate files from the tmp directory to avoid wasting space
	# NOTE: Review this periodically to make sure nothing's slipping through.
        # 
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes_plinked* 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.recode.vcf 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.matched.recode.vcf 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_prefiltered.recode_level{5}.vcf 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.log 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.log 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.gwas.config 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.eqtl.config 2> /dev/null".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)

if __name__ == "__main__":
	main()
