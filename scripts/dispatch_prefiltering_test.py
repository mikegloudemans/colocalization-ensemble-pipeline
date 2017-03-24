#!/usr/bin/python
import subprocess
import pandas as pd
import sys
import operator
import time
import random

gwas_threshold = 1e-7
window = 1000000

eqtl_threshold = "1e-05"
adaptive_threshold = "False"

gwas_file = sys.argv[1]
gwas_suffix = gwas_file.split("/")[-1].replace(".", "_")

threshold_directory = '{0:.1e}_by_{1:.1e}'.format(float(eqtl_threshold), gwas_threshold)
# Directory to be added to output path if not using adaptive thresholds
if not adaptive_threshold == "True":
        adapt_dir = "fixed_threshold/"
else:
        adapt_dir = ""

base_output_dir = "/users/mgloud/projects/brain_gwas/output/{0}/{1}{2}".format(gwas_suffix, adapt_dir, threshold_directory)

# Remove previous colocalization file
subprocess.call("rm {0}/{1}_clpp_status_prefiltering_test.txt".format(base_output_dir, gwas_suffix.replace(".", "_")), shell=True)

# Load in the data
gwas_table = pd.read_csv(gwas_file, sep="\t")

print gwas_table.head()

subset = gwas_table[['chr', 'snp_pos', 'pvalue']]

# For now, throw out all chrX results.
# I'm pretty sure the GTEx eQTL lists don't include chrX
subset = subset[subset['chr'] != 'chrX']
all_snps = [tuple(x) for x in subset.values]

all_snps = sorted(all_snps, key=operator.itemgetter(2))

snps_to_test = []
for snp in all_snps:
	if snp[2] >= gwas_threshold:
		break
	skip = False
	for kept_snp in snps_to_test:
		if kept_snp[0] == snp[0] and abs(kept_snp[1] - snp[1]) < window:
			skip = True
			break
	if not skip:
		snps_to_test.append(snp)

print sorted(snps_to_test)
print len(snps_to_test)

random.seed(2)
random.shuffle(snps_to_test)
snps_to_test = snps_to_test[:10]

for snp in snps_to_test:

	subprocess.check_call("python test_ecaviar_prefiltering.py {0} {1} {2} {3} {4} {5} {6} &".format(snp[0], snp[1], gwas_file, gwas_threshold, eqtl_threshold, adaptive_threshold, True), shell=True)
	subprocess.check_call("python test_ecaviar_prefiltering.py {0} {1} {2} {3} {4} {5} {6} &".format(snp[0], snp[1], gwas_file, gwas_threshold, eqtl_threshold, adaptive_threshold, False), shell=True)
	# Limit to 5 jobs at a time
	while int(subprocess.check_output('ps -ef | grep mgloud | grep "python test_ecaviar" | grep -v grep | grep {0} | wc -l'.format(gwas_file), shell=True)) > 5:
		time.sleep(1)
