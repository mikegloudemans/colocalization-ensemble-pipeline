#!/usr/bin/python
import subprocess
import pandas as pd
import sys
import operator
import time
import argparse
import os

def main():
	# Get command line input arguments
	args = parse_args()

	base_output_dir = get_output_dir(args)

	# Remove previous colocalization file
	subprocess.call("rm {0}/*clpp_status.txt".format(base_output_dir), shell=True)

	# Load in the data, get p-values for each SNP
	gwas_table = pd.read_csv(args.gwas_file, sep="\t")
	subset = gwas_table[['chr', 'snp_pos', 'pvalue']]

	# For now, throw out all chrX results.
	# I'm pretty sure the GTEx eQTL lists don't include chrX
	subset = subset[subset['chr'] != 'chrX']
	all_snps = [tuple(x) for x in subset.values]
	all_snps = sorted(all_snps, key=operator.itemgetter(2))

	snps_to_test = []
	for snp in all_snps:
		if snp[2] >= args.gwas_threshold:
			break

		# For now, ignore a SNP if it's in the MHC region, because there's way too
		# much signal there to make sense of it right now.
		if (str(snp[0]) == "6" or snp[0] == "chr6") and snp[1] > 25000000 and snp[1] < 35000000:
			continue

		skip = False
		for kept_snp in snps_to_test:
			if kept_snp[0] == snp[0] and abs(kept_snp[1] - snp[1]) < args.window:
				skip = True
				break
		if not skip:
			snps_to_test.append(snp)

	for snp in snps_to_test:
		print "python find_colocalization.py {0} {1} {2} {3} &".format(snp[0], snp[1], args.gwas_file, args.gwas_threshold)
		subprocess.check_call("python find_colocalization.py {0} {1} {2} {3} &".format(snp[0], snp[1], args.gwas_file, args.gwas_threshold), shell=True)
		# Limit to N jobs running on this GWAS file; user can specify desired number
		while int(subprocess.check_output('ps -ef | grep mgloud | grep "python find_coloc" | grep -v grep | grep {0} | wc -l'.format(args.gwas_file), shell=True)) >= args.max_jobs:
			time.sleep(1)

### Function parse_args
# Parses command line arguments
#
# Input: none
# Output: namespace containing command line arguments

def parse_args():
	# Parse input
	parser = argparse.ArgumentParser(description='Run colocalization pipeline for GWAS sites passing a given threshold')
	parser.add_argument("--gwas-threshold", dest="gwas_threshold", action="store", help="p-value cutoff for GWAS signficance")
	parser.add_argument("--gwas-file", dest="gwas_file", action = "store", help="file with GWAS summary statistics")
	parser.add_argument("--window", dest="window", action = "store", help="test for colocalization at SNPs within this distance from a GWAS hit", default=1000000)
	parser.add_argument("--max-jobs", dest="max_jobs", action = "store", help="max number of colocalization jobs to run at a time", default=5)
	args = parser.parse_args()

	# Validate input	
	try:
		args.gwas_threshold = float(args.gwas_threshold)
	except:
		print "Invalid threshold input. Must be a floating-point number."
		sys.exit()		
	try:
		args.window = int(args.window)
	except:
		print "Invalid window input. Must be an integer."
		sys.exit()
	try:
		args.max_jobs = int(args.max_jobs)
	except:
		print "Max jobs must be an interger."
		sys.exit()
		
	try:
		if not os.path.exists(args.gwas_file):
			print "Input file does not exist."
			sys.exit()
	except:
		print "No input file was specified."
		sys.exit()

	return args	

### Function get_output_dir
# Determine directory where we'll store outputs of the analysis
#
# Input: parsed command line arguments
# Output: directory for output files

def get_output_dir(args):
	gwas_suffix = args.gwas_file.split("/")[-1].replace(".", "_")
	threshold_directory = '{0:.1e}'.format(args.gwas_threshold)
	base_output_dir = "/users/mgloud/projects/brain_gwas/output/{0}/{1}".format(gwas_suffix, threshold_directory)

	return base_output_dir

if __name__ == "__main__":
	main()
