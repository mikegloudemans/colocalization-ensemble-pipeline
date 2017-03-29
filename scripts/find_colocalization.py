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

# TODO: Add ArgumentParser to parse all of this input directly, in a separate function
window = 1000000

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

	# Keep track of which GWAS-eQTL pairs have colocalized in which tissues
	coloc_status = {}
	all_tissues = []

	# Create base directory for output files, in case not already created
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/{0}".format(gwas_suffix), shell=True)
	subprocess.call("mkdir {0}".format(base_output_dir), shell=True)

	for tissue in os.listdir("/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_data/MatrixEQTL/allCisSNPGenePairs/"):

		tissue_prefix = "_".join(tissue.split(".")[0].split("_")[:-1])
		all_tissues.append(tissue_prefix)
		
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

		# Make directory for storing plots to visualize associations
		subprocess.call("mkdir {0}/gwas_by_eqtl_scatterplots".format(base_output_dir), shell=True)
		subprocess.call("mkdir {0}/gwas_by_eqtl_scatterplots/{1}_{2}".format(base_output_dir, gwas_chrom, gwas_pos), shell=True)
		subprocess.call("mkdir {0}/snp_overlaps".format(base_output_dir), shell=True)
		subprocess.call("mkdir {0}/snp_overlaps/{1}_{2}".format(base_output_dir, gwas_chrom, gwas_pos), shell=True)

		# TODO: This section is designed to run fastQTL to find conditional eQTLs, but currently
		# we've set the maximum conditional level to 0 so that it doesn't find conditional eQTLs.
		# Decide whether to fix this or just eliminate it if not necessary.
		
		for gene in genes:
			conditional_eqtls = eqtls[eqtls['gene'] == gene]
			covariates = []
			
			# TODO: Determine a better cutoff for when we'd consider doing
			# conditional eQTLs. I'm not sure if it's ever worth it though,
			# it kind of just overcomplicates things and we haven't found
			# anything from it yet.
			max_level = 0
			current_level = 0

			while current_level <= max_level:
				if current_level > 0:
					if not run_fastqtl(gene, covariates, tissue_prefix, gwas_suffix, gwas_chrom, gwas_pos, current_level):
						break

					with open("/users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_awked.txt".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, current_level, gwas_pos-window, gwas_pos+window), "w") as w:
						w.write("chr\tsnp_pos\tref\talt\tgenome\tgene\tbeta\tt-stat\tpvalue\n")
					subprocess.check_call('''sed s/_/" "/g /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt | awk '{{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $1 "\t" $10 "\t" $9 "\t" $8}}' >> /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_awked.txt'''.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, current_level, gwas_pos-window, gwas_pos+window), shell=True)

					# Load conditional eqtls from fastQTL results
					conditional_eqtls = \
						pd.read_csv('''/users/mgloud/projects/brain_gwas/tmp/fastQTL/''' + \
						'''{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_awked.txt'''.format(gwas_suffix, \
						gwas_chrom, gwas_pos, tissue_prefix, gene, current_level, gwas_pos-window, gwas_pos+window), sep="\t")

				if gwas_eqtl_colocalization(gwas_chrom, gwas_pos, gene, conditional_eqtls, \
					gwas_table, gwas_suffix, tissue_prefix, current_level) == False:
					break

				top_index = conditional_eqtls['pvalue'].tolist().index(min(conditional_eqtls['pvalue']))
				top = conditional_eqtls.iloc[top_index]
				
				covariates.append((str(gwas_chrom), str(top['snp_pos']), top['ref'], top['alt']))

				current_level += 1

				# For now, only go past conditional level 0 if it's a brain tissue.
				if "Brain" not in tissue:
					break
		
### Function gwas_eqtl_colocalization
# Test whether a gwas SNP and an eQTL colocalize. Currently will
# always return True unless it's impossible to perform the test.
# In future, might return True/False to signal whether we should
# run fastQTL again conditioned on the top remaining SNP.
#
# Input: gwas location, gene to test for eQTL colocalization, set of eqtls
# Output: True if test ran, False if test did not occur.
def gwas_eqtl_colocalization(gwas_chrom, gwas_pos, gene, eqtls, gwas_table, gwas_suffix, \
	tissue_prefix, current_level):

	# Get the set of genes we're interested in
	eqtl_subset = eqtls[eqtls['gene'] == gene]
	
	# TODO: View output logs to figure out why this is happening.
	# It shouldn't be happening, but apparently it is.
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

	# TODO: Modify plotting code so that eCAVIAR only plots results with colocalization, instead of results
	# passing a certain threshold. For now though, don't plot anything. (Could also delegate this 
	# to an entirely different Python script.

	#subprocess.call("mkdir {0}/gwas_by_eqtl_scatterplots/{1}_{2}".format(base_output_dir, gwas_chrom, gwas_pos), shell=True)
	#subprocess.call("mkdir {0}/gwas_by_eqtl_scatterplots/{1}_{2}/{3}".format(base_output_dir, gwas_chrom, gwas_pos, tissue_prefix), shell=True)

	#plt.figure(figsize=(10,10))
	#plt.scatter([-1 * math.log10(p) for p in not_passing['pvalue_x']], [-1 * math.log10(p) for p in not_passing['pvalue_y']])
	#plt.scatter([-1 * math.log10(p) for p in passing['pvalue_x']], [-1 * math.log10(p) for p in passing['pvalue_y']], color="red")
	#color_code = [int(sp)-gwas_pos for sp in combined['snp_pos']]
	#color_code = [min(abs(cc), 100000) * ((cc > 0)*2 - 1) for cc in color_code]
	#plt.scatter([-1 * math.log10(p) for p in combined['pvalue_x']], [-1 * math.log10(p) for p in combined['pvalue_y']], c=color_code, cmap=plt.cm.magma, edgecolor='', s=50)
	#plt.scatter([-1 * math.log10(p) for p in combined['pvalue_x']], [-1 * math.log10(p) for p in combined['pvalue_y']], c=combined['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)

	#if max([-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']]) < 20:
	#	plt.axis([0, 20, 0, 20])
	#else:
	#	plt.axis([0, max([20] + [-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']]), 0, max([20] + [-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']])])
	#plt.plot((-1*math.log10(gwas_threshold), -1*math.log10(gwas_threshold)), (0, max([20] + [-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']])), 'r--')
	#plt.plot((0, max([20] + [-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']])) ,(-1*math.log10(pval_threshold), -1*math.log10(pval_threshold)), 'r--')
	#plt.xlabel('GWAS -log p-value', fontsize=16)
	#plt.ylabel('eQTL -log p-value', fontsize=16)
	#plt.savefig("{0}/gwas_by_eqtl_scatterplots/{1}_{2}/{3}/{4}_level{5}.png".format(base_output_dir, gwas_chrom, gwas_pos, tissue_prefix, gene, current_level))
	#plt.gcf().clear()


	# Also create a LocusZoom-style plot showing the GWAS and eQTL signals next to one another.
	#plt.figure(figsize=(20,10))
	#plt.subplot(211)
	#plt.scatter(combined['snp_pos'], [-1 * math.log10(p) for p in combined['pvalue_x']], c=combined['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)
	#plt.plot((gwas_pos, gwas_pos), (-2, max([-1 * math.log10(p) for p in combined['pvalue_x']])), 'k--')
	#plt.ylabel('GWAS -log p-value', fontsize=16)
	#plt.subplot(212)
	#plt.scatter(combined['snp_pos'], [-1 * math.log10(p) for p in combined['pvalue_y']], c=combined['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)
	#plt.plot((gwas_pos, gwas_pos), (-2, max([-1 * math.log10(p) for p in combined['pvalue_y']])), 'k--')
	#plt.ylabel('eQTL -log p-value', fontsize=16)
	#plt.xlabel('Position', fontsize=16)
	#plt.savefig("{0}/gwas_by_eqtl_scatterplots/{1}_{2}/{3}/{4}_level{5}_lz.png".format(base_output_dir, gwas_chrom, gwas_pos, tissue_prefix, gene, current_level))
	#plt.gcf().clear()

	# Run eCAVIAR on SNPs.
	
	
	# This should make eCAVIAR run substantially faster. 
	# This is the newest version of subsetting.
	# Keep the union between the best 100 GWAS SNPs and best 100 eQTL SNPs,
	# since none of the others will have that much effect on the eCAVIAR score anyway.

	# TODO: Modify the analysis below to use FINEMAP instead of subsetting.
	# TODO: Run this test at a bunch of sample sites with FINEMAP and with eCAVIAR and
	# compare results to see if they're similar. If not, investigate why.
	combined = combined.sort_values(by='pvalue_x')
	ecaviar_x = combined.head(100)
	combined = combined.sort_values(by='pvalue_y')
	ecaviar_y = combined.head(100)

	ecaviar_set = pd.concat([ecaviar_x, ecaviar_y])
	ecaviar_set = ecaviar_set.drop_duplicates()
	run_ecaviar(ecaviar_set, gwas_chrom, gwas_pos, tissue_prefix, gene, current_level)

	
	return True


# Function: run_fastqtl	
# Returns True if fastqtl has been successfully run at the locus/gene pair.
# Returns False if we need to skip this locus/gene pair, probably due
# to lack of expression data. For now, we're not even using this function.
#
'''
def run_fastqtl(gene, covariates, tissue_prefix, gwas_suffix, gwas_chrom, gwas_pos, conditional_level, replace_existing = False):
	
	# First check output directory to see if fastQTL has already been run for this locus.
	# If so, then just skip this whole function.
	if not replace_existing:
		if os.path.exists("/users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)):
			return True

	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}".format(gwas_suffix.replace(".", "_")), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix), shell=True)

	# Get genotypes at covariates. For now, I'm going to use genotype likelihood estimates.
	# If this doesn't work I may need to figure out a different way to fill in the missing
	# values from the VCF.
	covariate_statuses = []
	for i in range(len(covariates)):
		covariate_statuses.append({})

	print covariates
	with open ("/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles/snpMatrices/{0}_Analysis.snps.txt".format(tissue_prefix)) as f:
		header = f.readline().strip().split()
		# NOTE/WARNING: If two possible variants exist at the same position, this approach will just pick the second one. 
		# Be wary of this; should fix this later to verify we're looking at the exact same variants as used in the 
		# original eQTL test.
		for line in f:
			data = line.strip().split()
			chr = data[0].split("_")[0]
			pos = data[0].split("_")[1]
			ref = data[0].split("_")[2]
			alt = data[0].split("_")[3]
			if (chr, pos, ref, alt) in covariates:
				ind = covariates.index((chr, pos, ref, alt))
				for i in range(1, len(data)):
					covariate_statuses[ind][header[i]] = data[i]
	print covariate_statuses

	# Create covariates file containing status for SNPs of interest
	with open("/users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_covariates_level{5}.txt".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), "w") as w:
		with open("/users/joed3/GTExCisEqtls/data/V6P/v6p_All_Tissues_covariates_FOR_QC_ONLY/{0}_Analysis.covariates.txt".format(tissue_prefix)) as f:
			header = f.readline()
			w.write(header)
			header = header.strip().split()

			# Rewrite the original file
			for line in f:
				if line.strip() != "":
					w.write(line)

			# Add lines at the end containing statuses of covariate SNPs.
			for i in range(len(covariates)):
				w.write("{0}_{1}\t".format(covariates[i][0], covariates[i][1]))
				for j in range(1,len(header)):
					w.write(covariate_statuses[i][header[j]] + "\t")
				w.write("\n")
				

	# NOTE: for the VCF subsetting section, we'll avoid subsetting for now if possible, and only do
	# that if it turns out things are really slow.

	# Subset VCF down to region of interest?

	# Subset gene expression data down to region of interest

	with gzip.open("/users/joed3/GTExCisEqtls/data/subsampling/{0}/{0}.normalized.expr.bed.gz".format(tissue_prefix), "rb") as f:
		with open("/users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_expression.bed".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene), "w") as w:
			w.write("#"+f.readline())

	# This isn't very elegant, but it's what we're doing for now. Run grep to see if the gene appears in the normalized expression list.
	# If not, then just skip the gene for now.
	# TODO: Fix this so that we're testing a consistent set of genes on all conditional levels (maybe throw away genes on the top level if
	# not found in the expression list).
	try:
		subprocess.check_call("zcat /users/joed3/GTExCisEqtls/data/subsampling/{3}/{3}.normalized.expr.bed.gz | grep {4} >> /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_expression.bed".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene), shell=True)
	except subprocess.CalledProcessError:
		return False

	# Check to make sure gene is within window of SNPs for fastQTL run. If not, we need to skip it.
        with open("/users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_expression.bed".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene)) as f:
                f.readline()
                loc = f.readline().strip().split()[1:3]

                if (loc[0] > gwas_pos + window and loc[1] > gwas_pos + window) or (loc[0] < gwas_pos - window and loc[1] < gwas_pos - window):
                        return False


	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_expression.bed.gz".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene), shell=True)
	subprocess.check_call("bgzip /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_expression.bed".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene), shell=True)
	subprocess.check_call("tabix -p bed /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_expression.bed.gz".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene), shell=True)

	# Run fastQTL
	subprocess.check_call("/users/zappala/software/fastqtl/bin/fastQTL \
				--vcf /users/joed3/GTExCisEqtls/data/subsampling/VCF/gtex.v6.allchr.impute.info04.maf01.hwep1e6.constrvarids.vcf.gz \
				--bed /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_expression.bed.gz \
				--window 1e6 \
				--maf-threshold 0.01 \
				--ma-sample-threshold 10 \
				--cov /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_covariates_level{5}.txt \
				--out /users/mgloud/projects/brain_gwas/tmp/fastQTL/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt \
				--region {1}:{6}-{7} \
				".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level, gwas_pos-window, gwas_pos+window), shell=True)
	return True
'''	

#def run_ecaviar(eqtl_chrom, eqtl_pos, gwas_pos, gene_list, gene, eqtls, gwas_table):
def run_ecaviar(combined, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level):

	# Make required directories for eCAVIAR analysis
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/".format(gwas_suffix), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}".format(gwas_suffix, gwas_chrom, gwas_pos), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/plink/{0}/".format(gwas_suffix), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}".format(gwas_suffix, gwas_chrom, gwas_pos), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/".format(gwas_suffix), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}".format(gwas_suffix, gwas_chrom, gwas_pos), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), shell=True)
	
	# This step is critical. If we don't call it, then eCAVIAR will append results to the 
	# previously generated list, which will probably screw everything up.
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
	# TODO: For the sake of a speedy analysis, maybe require MAF above 0.05, since hard to find eQTLs otherwise.
	# Only if necessary though.
	command = 'vcftools --vcf /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_prefiltered.recode_level{5}.vcf --positions /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt --recode --recode-INFO-all --out /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes'.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)
        subprocess.check_call(command, shell=True)

	# TODO: Consider whether we should be using GTEx VCFs to compute LD instead of
	# using 1000 Genomes VCFs. Problem there is it's tough to do the
	# computations, since many individuals are duplicated in these VCFs. Would have to
	# get consensus genotypes for each person based on samples across all tissues, or something like that.

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
	subprocess.check_call("sed s/nan/0/g /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.ld > /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)

	# Print z-scores to input files for eCAVIAR
	combined['ZSCORE_eqtl'] = combined['t-stat']
	# Figure out whether GWAS scores are in odds ratio or beta-se format
	if 'or' in combined:
		combined['ZSCORE_gwas'] = (combined['or']-1) / combined['se']
	else:
		combined['ZSCORE_gwas'] = (combined['beta_x']) / combined['se']
	with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), "w") as w:
		snps = combined[['snp_pos', 'ZSCORE_eqtl']]
		snps.to_csv(w, index=False, header=False, sep="\t")

 	with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), "w") as w:
		snps = combined[['snp_pos', 'ZSCORE_gwas']]	
		snps.to_csv(w, index=False, header=False, sep="\t")

	# Run eCAVIAR
	# NOTE: Redirecting output to /dev/null because eCAVIAR prints the entire
	# LD matrix, which is just bad when you have thousands of loci. Fix this if
	# you need to see eCAVIAR output.
	command = '''/srv/persistent/bliu2/tools/caviar/CAVIAR-C++/eCAVIAR \
		-o /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_ecaviar_results \
		-l /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld \
		 -l /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld \
		 -z /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z \
		  -z /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z \
		-c 1 \
		-r 0.95 > /dev/null'''.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)
	subprocess.check_call(command, shell=True)

	# Parse eCAVIAR results to compute CLPP score
	command = '''awk '{{sum += $2}} END {{print sum}}' ../tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_ecaviar_results_col'''.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)
	clpp = float(subprocess.check_output(command, shell=True))

	# Add results to the desired file
	with open("{0}/{1}_clpp_status.txt".format(base_output_dir, gwas_suffix.replace(".", "_")), "a") as a:
		a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level, snps.shape[0], clpp))

	# Remove temporary intermediate files to save space
	purge_tmp_files(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level)

### Function purge_tmp_files
# Remove temporary files created during this run of eCAVIAR,
# to free up space.
def purge_tmp_files(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level):
	
	# Why not just purge everything from the entire directory of the GWAS that we're working on?
	# Answer: it's possible that other jobs are still using those files.
	# Only purge the specific files created in this run.

	# Remove intermediate files from the tmp directory to avoid wasting space
	# NOTE: Review this periodically to make sure nothing's slipping through.
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.txt".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.fixed.ld".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes_plinked*".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.recode.vcf".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.matched.recode.vcf".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}/{4}_prefiltered.recode_level{5}.vcf".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)
	subprocess.call("rm /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_1Kgenomes.log".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix, gene, conditional_level), shell=True)


if __name__ == "__main__":
	main()
