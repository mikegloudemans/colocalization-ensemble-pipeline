# Author: Mike Gloudemans
#
# preprocess.py
#
# Tools for loading summary statistics.
#

import subprocess
import pandas as pd
import operator

import sys
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO


# Input: gwas file location, threshold of significance,
#   minimum distance between two selected SNPs.
# Output: A list of significant SNPs to test.
def select_test_snps(gwas_file, gwas_threshold, window=1000000):

    print("Selecting GWAS hits from {0}".format(gwas_file))

    # Load in the data, get p-values for each SNP
    gwas_table = pd.read_csv(gwas_file, sep="\t")
    subset = gwas_table[['chr', 'snp_pos', 'pvalue']]

    # For now, throw out all chrX results.
    subset = subset[subset['chr'] != 'chrX']
    all_snps = [tuple(x) for x in subset.values]
    all_snps = sorted(all_snps, key=operator.itemgetter(2))

    # Go through the list of SNPs in order, adding the ones
    # passing our criteria.
    snps_to_test = []
    for snp in all_snps:

        # See if we're done yet        
        if snp[2] >= gwas_threshold:
                break

        # For now, ignore a SNP if it's in the MHC region -- this
        # would require alternative methods.
        if (str(snp[0]) == "6" or snp[0] == "chr6") and snp[1] > 25000000 and snp[1] < 35000000:
                continue

        # Before adding a SNP, make sure it's not right next
        # to another SNP that we've already selected.
        skip = False
        for kept_snp in snps_to_test:
                if kept_snp[0] == snp[0] and abs(kept_snp[1] - snp[1]) < window:
                        skip = True
                        break
        if not skip:
                snps_to_test.append(snp)

    return snps_to_test


# Load summary statistics for GWAS
def get_gwas_data(gwas_file, snp, window=500000):

    gwas_chrom = snp[0]
    gwas_pos = snp[1]

    # Subset GWAS list to SNPs near the GWAS position
    gwas_table = pd.read_csv(gwas_file, sep="\t")
    gwas_table = gwas_table[(gwas_table['snp_pos'] > gwas_pos - window) & (gwas_table['snp_pos'] < gwas_pos + window)]
    gwas_table = gwas_table[(gwas_table['chr'] == gwas_chrom) | (gwas_table['chr'] == 'chr{0}'.format(gwas_chrom))]

    return gwas_table

# Load summary statistics for eQTL
def get_eqtl_data(eqtl_file, snp, window=500000):

    gwas_chrom = snp[0][3:]
    gwas_pos = snp[1]

    # Get eQTL data using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 1".format(eqtl_file), shell=True)
    raw_eqtls = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(eqtl_file, \
            gwas_chrom, gwas_pos - window, gwas_pos + window), shell=True)
    eqtls = pd.read_csv(StringIO(header + raw_eqtls), sep="\t")

    return eqtls

# Input: GWAS pandas dataframe, eQTL pandas dataframe, gene name as a string,
#   GWAS SNP as a tuple.
# Returns: a combined table of summary statistics, or None if we need to skip
#   the site due to insufficient data.
def combine_summary_statistics(gwas_data, eqtl_data, gene, snp):

    gwas_pos = snp[1]

    # Filter SNPs down to the gene of interest.
    eqtl_subset = eqtl_data[eqtl_data['gene'] == gene]
	
    # Sometimes the GWAS SNP is outside of the range of eQTLs tested for a certain
    # gene, or on the outside fringe of the range. If this is the case, then skip it.
    # NOTE: Modify the 50000 cutoff if it doesn't seem like it's giving enough room for LD decay to fall off.

    if gwas_pos > max(eqtl_subset['snp_pos']) + 50000 or gwas_pos < min(eqtl_subset['snp_pos'] - 50000):
            return None
    
    # Join the list of eQTL SNPs with the list of GWAS SNPs
    combined = pd.merge(gwas_data, eqtl_subset, on="snp_pos")

    # Check to make sure there are SNPs remaining; if not, just move on
    # to next gene.
    if combined.shape[0] == 0: 
            return None

    return combined




