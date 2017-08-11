# Author: Mike Gloudemans
#
# preprocess.py
#
# Tools for loading summary statistics.
#

import subprocess
import pandas as pd
import operator
import SNP
from scipy import stats


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

    all_snps = [tuple(x) for x in subset.values]
    all_snps = sorted(all_snps, key=operator.itemgetter(2))

    # For now, include only autosomal SNPs.
    filtered = []
    for s in all_snps:
        if "chr" in str(s[0]):
            try:
                filtered.append((int(s[0][3:]), s[1], s[2]))
            except:
                pass
        else:
            try:
                filtered.append((int(s[0]), s[1], s[2]))
            except:
                pass

    all_snps = [SNP.SNP(x) for x in filtered]

    # Go through the list of SNPs in order, adding the ones
    # passing our criteria.
    snps_to_test = []
    for snp in all_snps:

        # See if we're done yet        
        if snp.pval >= gwas_threshold:
                break

        # For now, ignore a SNP if it's in the MHC region -- this
        # would require alternative methods.
        if (snp.chrom == 6) and snp.pos > 25000000 and snp.pos < 35000000:
                continue

        # Before adding a SNP, make sure it's not right next
        # to another SNP that we've already selected.
        skip = False
        for kept_snp in snps_to_test:
                if kept_snp.chrom == snp.chrom and abs(kept_snp.pos - snp.pos) < window:
                        skip = True
                        break
        if not skip:
                snps_to_test.append(snp)

    
    print("Testing {0} SNPs.".format(len(snps_to_test)))

    return snps_to_test


# Load summary statistics for GWAS
def get_gwas_data(gwas_file, snp, window=500000):

    # Subset GWAS list to SNPs near the GWAS position
    gwas_table = pd.read_csv(gwas_file, sep="\t")
    gwas_table['snp_pos'] = gwas_table['snp_pos'].astype(int)
    gwas_table = gwas_table[(gwas_table['snp_pos'] > snp.pos - window) & (gwas_table['snp_pos'] < snp.pos + window)]
    gwas_table = gwas_table[(gwas_table['chr'] == snp.chrom) | (gwas_table['chr'] == 'chr{0}'.format(snp.chrom))]


    # Figure out whether GWAS scores are in odds ratio or beta-se format
    if 'or' in gwas_table:
        gwas_table['ZSCORE'] = (gwas_table['or']-1) / gwas_table['se']
    elif 'beta_x' in gwas_table:
        gwas_table['ZSCORE'] = (gwas_table['beta_x']) / gwas_table['se']
    elif 'pvalue' in gwas_table and "direction" in gwas_table:
        # TODO: Test this
        gwas_table['ZSCORE'] = stats.norm.isf(gwas_table["pvalue"] / 2) * (2*(gwas_table["direction"] == "+")-1)
    else:
        return None


    return gwas_table

# Load summary statistics for eQTL
def get_eqtl_data(eqtl_file, snp, window=500000):

    # Get eQTL data using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 1".format(eqtl_file), shell=True)
    raw_eqtls = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(eqtl_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True)
    eqtls = pd.read_csv(StringIO(header + raw_eqtls), sep="\t")
    eqtls['snp_pos'] = eqtls['snp_pos'].astype(int)

    if 't-stat' in eqtls:
        eqtls['ZSCORE'] = eqtls['t-stat']
    elif "chisq" in eqtls:
        # TODO: Test this because I'm not sure yet if it's correct.
        eqtls['pvalue'] = 1-stats.chi2.cdf(eqtls["chisq"],1)
        eqtls['ZSCORE'] = stats.norm.isf(eqtls['pvalue']/2) * (2 * (eqtls["effect_size"] > 0) - 1) # Might not be correct: all have effect size > 0 right now...
    else:
        return None



    return eqtls

# Input: GWAS pandas dataframe, eQTL pandas dataframe, gene name as a string,
#   GWAS SNP as a tuple.
# Returns: a combined table of summary statistics, or None if we need to skip
#   the site due to insufficient data.
def combine_summary_statistics(gwas_data, eqtl_data, gene, snp):

    # Filter SNPs down to the gene of interest.
    eqtl_subset = eqtl_data[eqtl_data['gene'] == gene]
	
    # Sometimes the GWAS SNP is outside of the range of eQTLs tested for a certain
    # gene, or on the outside fringe of the range. If this is the case, then skip it.
    # NOTE: Modify the 50000 cutoff if it doesn't seem like it's giving enough room for LD decay to fall off.

    if snp.pos > max(eqtl_subset['snp_pos']) + 50000 or snp.pos < min(eqtl_subset['snp_pos'] - 50000):
            return None
 
    # Skip it if there's nothing left
    if gwas_data.shape[0] == 0:
            return None

    # For now, filter out sites where p-values are too extreme.
    # TODO: Modify the way we handle these so we can still use extremely
    # significant sites.

    if min(eqtl_subset['pvalue']) < 1e-300 or min(gwas_data['pvalue']) < 1e-300:
            return None

    # TODO TODO: At this step filter out multi-allelic sites, any variants
    # whose same position appears twice or more. (can steal the code 
    # that is currently used in the FINEMAP pipeline). This is important 
    # because otherwise at this level some sites with 2 SNPs might be merged into
    # 2x2 SNPs or something like that

    # Join the list of eQTL SNPs with the list of GWAS SNPs
    combined = pd.merge(gwas_data, eqtl_subset, on="snp_pos", suffixes=("_gwas", "_eqtl"))

    # Check to make sure there are SNPs remaining; if not, just move on
    # to next gene.
    if combined.shape[0] == 0: 
            return None

    return combined

# Use vcftools to get MAFs for variants in a dataset
# NOTE: This function would throw an error if there's a variant
# in our lists that's not also in 1K genomes. For now though, let's
# just cross that bridge when we get to it. It will probably come
# down to a simple fix of just throwing away the variants that 
# don't appear in 1000 genomes, since they're low-powered anyway.
def add_maf(data):

    # TODO: Under construction.

    # Write variants to a file, formatted for VCFtools
    pd.to_csv("some_temp_file")

    # Run VCFtools to 

    # Parse results from VCFtools, and add them to our data

    # Remove both temp files

    return data


