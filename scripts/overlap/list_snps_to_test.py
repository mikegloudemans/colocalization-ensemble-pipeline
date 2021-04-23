#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 6/7/2018
# Modified: 4/23/2021

import operator
import gzip
import sys
import subprocess

hits_file = sys.argv[1]
lookup_file = sys.argv[2]
output_file = sys.argv[3]
lookup_threshold = float(sys.argv[4])
lookup_window = int(sys.argv[5])

def main():

	print(f"Getting overlaps for {hits_file} and {lookup_file}.")

	with gzip.open(lookup_file) as f:
		header = f.readline().decode("utf-8").strip().split()
		header = [h.lower() for h in header]

	pval_index = header.index("pvalue")

	# Is there a column specifying the trait / gene / feature in the lookup target file?
	# If so, get the index of it
	if "trait" in header:
		lookup_trait_index = header.index("trait")
	elif "gene" in header:
		lookup_trait_index = header.index("gene")
	elif "feature" in header:
		lookup_trait_index = header.index("feature")
	else:
		lookup_trait_index = -1
  
	with open(output_file, "w") as w:

		w.write("chrom\tpos\tsource_pvalue\tsource_trait\tlookup_pvalue\tlookup_trait\n")

		with open(hits_file) as f:
			f.readline()
			for line in f:
				chrom, pos, source_pvalue, source_trait = line.strip().split()

				min_pos = int(pos)-lookup_window
				max_pos = int(pos)+lookup_window

				wide_matches = subprocess.run(f"tabix {lookup_file} chr{chrom}:{min_pos}-{max_pos}".split(), capture_output=True).stdout.decode('utf-8')
				wide_matches += subprocess.run(f"tabix {lookup_file} {chrom}:{min_pos}-{max_pos}".split(), capture_output=True).stdout.decode('utf-8')

				if wide_matches == "":
					continue
	
				# Sort by pval so we can be sure we get the most significant SNP at the locus first
				wide_matches = wide_matches.strip().split("\n")
				wide_matches = [wm.split("\t") for wm in wide_matches]

				def numerize_pval(x):
					new = x[:]
					new[pval_index] = float(new[pval_index])
					return(new)

				try:
					wide_matches = [numerize_pval(d) for d in wide_matches]
				except:
					print("Formatting error: could not convert p-value to float")
					return
				wide_matches = sorted(wide_matches, key=operator.itemgetter(pval_index))

				matched = set([])	 # Don't print repeats if there are multiple matches				
				for data in wide_matches:

					# If there's more than one trait in the lookup file,
					# figure out which trait this row is measuring
					if lookup_trait_index != -1:
						lookup_trait = data[lookup_trait_index]
					else:
						lookup_trait = "none"

					# Keep track of SNP passing out overlap cutoff in the lookup set

					lookup_pvalue = data[pval_index]

					if float(data[pval_index]) <= lookup_threshold and lookup_trait not in matched:
						w.write(f"{chrom}\t{pos}\t{source_pvalue}\t{source_trait}\t{lookup_pvalue}\t{lookup_trait}\n")
						matched.add(lookup_trait)

if __name__ == "__main__":
	main()
