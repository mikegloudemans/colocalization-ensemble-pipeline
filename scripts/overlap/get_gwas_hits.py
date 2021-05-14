#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 6/7/2018
# Modified: 4/14/2021

# Get all independent GWAS hits for the input file
# that pass a certain threshold. Output to a temporary
# file.

import sys
import operator
import gzip

gwas_file = sys.argv[1]
out_file = sys.argv[2]
gwas_threshold = float(sys.argv[3])
gwas_window = int(sys.argv[4])

def main():
	with open(out_file, "w") as w:
		w.write("chr\tsnp_pos\tpvalue\ttrait\n")

		hits = get_gwas_hits(gwas_file, gwas_threshold, gwas_window)
		for hit in hits:
			print(hit)
			w.write("\t".join([str(h) for h in hit])+"\n")
			

def get_gwas_hits(gwas_file, gwas_threshold, gwas_window):

	print(f"Getting GWAS hits for {gwas_file}")
	print(f"P < {gwas_threshold}")
	print(f"At least {gwas_window} BP apart")
	print()

	with gzip.open(gwas_file) as f:

		header = f.readline().strip().split()
		header = [str(h.decode().lower()) for h in header]


		if "trait" in header:
			source_trait_index = header.index("trait")
		elif "gene" in header:
			source_trait_index = header.index("gene")
		elif "feature" in header:
			source_trait_index = header.index("feature")
		else:
			source_trait_index = -1

		pval_index = header.index("pvalue")
		chr_index = header.index("chr")
		snp_pos_index = header.index("snp_pos")

		all_snps = []

		trait = "none"
		for line in f:
			data = line.decode("utf-8").strip().split("\t")
			if source_trait_index != -1:
				trait = data[source_trait_index]
			try:
				pvalue = float(data[pval_index])
			except:
				continue
			chr = data[chr_index]
			pos = int(data[snp_pos_index])
			if pvalue > gwas_threshold:
				continue

			all_snps.append((chr, pos, pvalue, trait))
			
	# For now, include only autosomal SNPs.
	filtered = []
	for s in all_snps:
		if "chr" in str(s[0]):
			filtered.append((s[0][3:], s[1], s[2], s[3]))
		else:
			filtered.append((s[0], s[1], s[2], s[3]))

	all_snps = sorted(filtered, key=operator.itemgetter(2)) 

	# Go through the list of SNPs in order, adding the ones
	# passing our criteria.
	snps_to_test = []
	for snp in all_snps:

		# For now, ignore a SNP if it's in the MHC region -- this
		# would require alternative methods.
		if (snp[0] == "6") and snp[1] > 25000000 and snp[1] < 35000000:
			continue

		# Before adding a SNP, make sure it's not right next
		# to another SNP that we've already selected.
		
		skip = False
		for kept_snp in snps_to_test:
			if kept_snp[0] == snp[0] and abs(kept_snp[1] - snp[1]) < gwas_window and kept_snp[3] == snp[3]:
				skip = True
				break
		if not skip:
			snps_to_test.append(snp)

	return snps_to_test

if __name__ == "__main__":
	main()
