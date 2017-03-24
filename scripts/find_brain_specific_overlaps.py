#!/usr/bin/python

import os
import sys

for snp in os.listdir("/users/mgloud/projects/brain_gwas/output/coloc_status/scz2_snp_results_txt/chunks"):
	with open("/users/mgloud/projects/brain_gwas/output/coloc_status/scz2_snp_results_txt/chunks/{0}".format(snp)) as f:
		header = f.readline().strip().split()
		for line in f:
			data = line.strip().split()
			snp = data[0]
			gene = data[1]
			brain_count = 0
			nonbrain_count = 0
			for i in range(2, len(data)):
				if data[i] == "1":
					if "Brain" in header[i]:
						brain_count += 1
					else:
						nonbrain_count += 1
			if brain_count > nonbrain_count:
				print snp, gene				
