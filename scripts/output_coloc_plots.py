#!/usr/bin/python
# Author: Mike Gloudemans
# Date: 11/16/2016

# Command-line input: 
# ./output_coloc_plots.py gwas_chrom gwas_pos gwas_file eqtl_gene
#
# Will create scatterplots in every tissue, for this pair
# of GWAS SNP and eQTL.
#
# Note: could probably be rolled back together with the find_colocalization
# script at some point.
#

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

gwas_chrom = int(sys.argv[1].replace('chr', ''))
gwas_pos = int(sys.argv[2])
gwas_file = sys.argv[3]
eqtl_gene = sys.argv[4]

pval_threshold = 1e-7

window = 1000000

def main():
	# Subset GWAS list to SNPs within 1MB of the GWAS position
	
	gwas_table = pd.read_csv(gwas_file, sep="\t")

	gwas_table = gwas_table[(gwas_table['snp_pos'] > gwas_pos - window) & (gwas_table['snp_pos'] < gwas_pos + window)]
	gwas_table = gwas_table[(gwas_table['chr'] == gwas_chrom) | (gwas_table['chr'] == 'chr{0}'.format(gwas_chrom))]
	
	gwas_prefix = "/".join(gwas_file.split("/")[:-1])
	gwas_suffix = gwas_file.split("/")[-1]

	# Subset eQTL list to SNPs with 1MB of the GWAS position.

	all_tissues = []
	for tissue in os.listdir("/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_data/MatrixEQTL/allCisSNPGenePairs/"):


		tissue_prefix = "_".join(tissue.split(".")[0].split("_")[:-1])
		all_tissues.append(tissue_prefix)


		# If eqtl file already exists, go straight to it. Otherwise, we need to do the time-consuming step
		# of creating it.
		try:
			eqtls = pd.read_csv("/users/mgloud/projects/brain_gwas/data/eqtls/{0}/{1}_{2}/{3}.txt".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), sep="\t")
	
		except:
			subprocess.call("mkdir /users/mgloud/projects/brain_gwas/data/eqtls/{0}".format(gwas_suffix), shell=True)
			subprocess.call("mkdir /users/mgloud/projects/brain_gwas/data/eqtls/{0}/{1}_{2}".format(gwas_suffix, gwas_chrom, gwas_pos), shell=True)
			# This is really gnarly. I'm sorry for anyone who has to read this
			subprocess.check_call('''echo "chr\tsnp_pos\tref\talt\tgenome\tgene\tbeta\tt-stat\tpvalue" > /users/mgloud/projects/brain_gwas/tmp/eqtl_subsets/{0}_{1}_{2}_header.txt'''.format(gwas_suffix, gwas_chrom, gwas_pos), shell=True)
			subprocess.check_call('''zcat /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_data/MatrixEQTL/allCisSNPGenePairs/{5} | sed 's/_/\t/g' | awk '{{if ($1 == {0} && $2 > {1} && $2 < {2}) print $0}}' > /users/mgloud/projects/brain_gwas/tmp/eqtl_subsets/{3}_{0}_{4}.tail.tmp'''.format(gwas_chrom, gwas_pos - window, gwas_pos + window, gwas_suffix, gwas_pos, tissue), shell=True)
			subprocess.check_call('''cat /users/mgloud/projects/brain_gwas/tmp/eqtl_subsets/{0}_{1}_{2}_header.txt /users/mgloud/projects/brain_gwas/tmp/eqtl_subsets/{0}_{1}_{2}.tail.tmp > /users/mgloud/projects/brain_gwas/data/eqtls/{0}/{1}_{2}/{3}.txt'''.format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), shell=True)
			subprocess.check_call('''rm /users/mgloud/projects/brain_gwas/tmp/eqtl_subsets/{0}_{1}_{2}_header.txt'''.format(gwas_suffix, gwas_chrom, gwas_pos), shell=True)
			subprocess.check_call('''rm /users/mgloud/projects/brain_gwas/tmp/eqtl_subsets/{0}_{1}_{2}.tail.tmp'''.format(gwas_suffix, gwas_chrom, gwas_pos), shell=True)

			# Read in relevant eQTLs
			eqtls = pd.read_csv("/users/mgloud/projects/brain_gwas/data/eqtls/{0}/{1}_{2}/{3}.txt".format(gwas_suffix, gwas_chrom, gwas_pos, tissue_prefix), sep="\t")

		# Get full list of all genes whose eQTLs we're testing at this site
		genes = set(eqtls['gene'])

		# Make directory for storing plots to visualize associations
		subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/gwas_by_eqtl_scatterplots/{0}".format(gwas_suffix.replace(".", "_")), shell=True)
		subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/gwas_by_eqtl_scatterplots/{0}/{1}_{2}".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos), shell=True)
		subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/snp_overlaps/{0}".format(gwas_suffix.replace(".", "_")), shell=True)
		subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/snp_overlaps/{0}/{1}_{2}".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos), shell=True)

		# TODO: Keep a list of all the SNP-gene pairs we've tested to get an idea of the scale of this.
		for gene in genes:
			gene_by_eqtl_comparison(gwas_chrom, gwas_pos,  gene, eqtls, gwas_table, gwas_suffix, coloc_status, tissue_prefix)

		subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/coloc_status/", shell=True)
		subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/coloc_status/{0}".format(gwas_suffix.replace(".", "_")), shell=True)
		subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/coloc_status/{0}/chunks".format(gwas_suffix.replace(".", "_")), shell=True)

	with open("/users/mgloud/projects/brain_gwas/output/coloc_status/{0}/chunks/{1}_{2}.txt".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos), "w") as w:
		w.write("snp\tgene\t{0}\n".format("\t".join(sorted(all_tissues))))
		for gene in coloc_status.keys():
			w.write("{0}_{1}\t{2}".format(gwas_chrom, gwas_pos, gene))
			for t in sorted(all_tissues):
				if t in coloc_status[gene]:
					w.write("\t1")
				else:
					w.write("\t0")
					
			w.write("\n")



def gene_by_eqtl_comparison(gwas_chrom, gwas_pos, gene, eqtls, gwas_table, gwas_suffix, coloc_status, tissue_prefix):

	if gene not in coloc_status:
		coloc_status[gene] = set([])

	# Get the set of genes we're interested in
	eqtl_subset = eqtls[eqtls['gene'] == gene]

	# Join the list of eQTL SNPs with the list of GWAS SNPs
	combined = pd.merge(gwas_table, eqtl_subset, on="snp_pos")

	# Check to make sure there are SNPs remaining; if not, just move on
	# to next gene.
	if combined.shape[0] == 0: 
		return

	# Subset down to SNPs that have reached reasonable significance in both datasets

	passing = combined[(combined['pvalue_x'] < pval_threshold) & (combined['pvalue_y'] < pval_threshold)]
	not_passing = combined[(combined['pvalue_x'] >= pval_threshold) | (combined['pvalue_y'] >= pval_threshold)]
	
	# If any SNPs coincide between the two datasets, make a note of it and save a data table with the overlapping snps,
	# as well as a scatterplot.
	# If not, then just move on.
	if passing.shape[0] == 0:
		return


	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/gwas_by_eqtl_scatterplots/{0}/{1}_{2}".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos), shell=True)
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/gwas_by_eqtl_scatterplots/{0}/{1}_{2}/{3}".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix), shell=True)

	plt.scatter([-1 * math.log10(p) for p in not_passing['pvalue_x']], [-1 * math.log10(p) for p in not_passing['pvalue_y']])
	plt.scatter([-1 * math.log10(p) for p in passing['pvalue_x']], [-1 * math.log10(p) for p in passing['pvalue_y']], color="red")
	if max([-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']]) < 30:
		plt.axis([0, 30, 0, 30])
	else:
		plt.axis([0, max([30] + [-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']]), 0, max([30] + [-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']])])
	plt.plot((-1*math.log10(pval_threshold), -1*math.log10(pval_threshold)), (0, max([30] + [-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']])), 'r--')
	plt.plot((0, max([30] + [-1 * math.log10(p) for p in combined['pvalue_x']] + [-1 * math.log10(p) for p in combined['pvalue_y']])) ,(-1*math.log10(pval_threshold), -1*math.log10(pval_threshold)), 'r--')
	plt.xlabel('GWAS log p-value', fontsize=16)
	plt.ylabel('eQTL log p-value', fontsize=16)
	plt.savefig("/users/mgloud/projects/brain_gwas/output/gwas_by_eqtl_scatterplots/{0}/{1}_{2}/{3}/{4}.png".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene))
	plt.gcf().clear()


	
	subprocess.call("mkdir /users/mgloud/projects/brain_gwas/output/snp_overlaps/{0}/{1}_{2}/{3}".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix), shell=True)
	with open("/users/mgloud/projects/brain_gwas/output/snp_overlaps/{0}/{1}_{2}/{3}/{4}.txt".format(gwas_suffix.replace(".", "_"), gwas_chrom, gwas_pos, tissue_prefix, gene), "w") as w:
		passing.to_csv(w, index=False, sep="\t")

	coloc_status[gene].add(tissue_prefix)
		



if __name__ == "__main__":
	main()
