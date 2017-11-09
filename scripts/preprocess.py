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
import math

import sys 
if sys.version_info[0] < 3: 
    from StringIO import StringIO 
else: 
    from io import StringIO


# Input: gwas file location, threshold of significance, minimum distance
# between two selected SNPs.  Output: A list of significant SNPs to test.
def select_test_snps(gwas_file, gwas_threshold, window=1000000):

    print("Selecting GWAS hits from {0}".format(gwas_file))

    stream = StringIO(subprocess.check_output("zcat {0}".format(gwas_file), shell=True))
    gwas_table = pd.read_csv(stream, sep="\t")
    subset = gwas_table[['chr', 'snp_pos', 'pvalue']]
    # TODO: Fix this line! Something is wrong with it I guess
    subset['pvalue'] = subset['pvalue'].astype(float)
    subset = subset[subset['pvalue'] < gwas_threshold]

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

    # TODO: Tabix this to make it faster, or at least
    # load the table just once for the locus. Right now
    # it's really slow.

    # Subset GWAS list to SNPs near the GWAS position
    '''gwas_table = pd.read_csv(gwas_file, sep="\t")
    gwas_table['snp_pos'] = gwas_table['snp_pos'].astype(int)
    gwas_table = gwas_table[(gwas_table['snp_pos'] > snp.pos - window) & (gwas_table['snp_pos'] < snp.pos + window)]
    gwas_table = gwas_table[(gwas_table['chr'] == snp.chrom) | (gwas_table['chr'] == 'chr{0}'.format(snp.chrom))]'''

    # Get GWAS data using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 1".format(gwas_file), shell=True)
    raw_gwas = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(gwas_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True)
    gwas_table = pd.read_csv(StringIO(header + raw_gwas), sep="\t")
    gwas_table['snp_pos'] = gwas_table['snp_pos'].astype(int)

    # Figure out whether GWAS scores are in odds ratio or beta-se format
    # NOTE: This section is likely to be error prone at the moment...be careful!
    # TODO: Eliminate options here and require a standard preprocessing format
    # to make things work more smoothly
    if 'or' in gwas_table:
        # TODO: Also verify the correctness of this math. Is it right?
        gwas_table['ZSCORE'] = (gwas_table['or']-1) / gwas_table['se']
    elif 'beta' in gwas_table:
        gwas_table['ZSCORE'] = (gwas_table['beta']) / gwas_table['se']
    elif 'pvalue' in gwas_table and "direction" in gwas_table:
        # This is the format for RPE.
        # Need to cap it at z-score of 40 for outrageous p-values (like with AMD / RPE stuff)        
        gwas_table['ZSCORE'] = pd.Series([min(x, 40) for x in stats.norm.isf(gwas_table["pvalue"] / 2)]) * (2*(gwas_table["direction"] == "+")-1)

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

    print eqtls.head()

    eqtls['snp_pos'] = eqtls['snp_pos'].astype(int)

    if 't-stat' in eqtls:
        eqtls['ZSCORE'] = eqtls['t-stat']
    elif "chisq" in eqtls:
        # TODO: Fix this to make it more explicitly clear that we're dealing with RASQUAL data
        # where effect size is given by allelic imbalance percentage pi
        # Use max function to protect against underflow in chi2 computation
        # TODO: Fixing z-score issues here is VERY important! This could cause serious issues downstream...
        # Do this even before starting on ASE stuff!
        # TODO TODO TODO
        eqtls['pvalue'] = [max(x, 1e-16) for x in 1-stats.chi2.cdf(eqtls["chisq"],1)]
        eqtls['ZSCORE'] = stats.norm.isf(eqtls['pvalue']/2) * (2 * (eqtls["pi"] > 0.5) - 1)
    else:
        return None



    return eqtls

# Input: GWAS pandas dataframe, eQTL pandas dataframe, gene name as a string,
#   GWAS SNP as a tuple.
# Returns: a combined table of summary statistics, or None if we need to skip
#   the site due to insufficient data.
def combine_summary_statistics(gwas_data, eqtl_data, gene, snp, settings, window=500000, unsafe=False):

    # TODO TODO TODO: Fix the SettingWithCopyWarning. It seems likely to be error-prone, according
    # to the manual!

    # Filter SNPs down to the gene of interest.
    eqtl_subset = eqtl_data[eqtl_data['gene'] == gene]

    # Sometimes the GWAS SNP is outside of the range of eQTLs tested for a certain
    # gene, or on the outside fringe of the range. If this is the case, then skip it.
    # NOTE: Modify the 50000 cutoff if it doesn't seem like it's giving enough room for LD decay to fall off.

    if snp.pos > max(eqtl_subset['snp_pos']) + 50000 or snp.pos < min(eqtl_subset['snp_pos'] - 50000):
            return "SNP outside range."
 
    # If not explicitly allowing them, remove pvalues with danger
    # of underflow.
    if min(eqtl_subset['pvalue']) < 1e-150:
        if unsafe:
            eqtl_subset['pvalue'] = eqtl_subset['pvalue'].apply(lambda x: max(x, 1e-150))
        else:
            return "eQTL pvalue underflow."

    if min(gwas_data['pvalue']) < 1e-150:
        if unsafe:
            gwas_data['pvalue'] = gwas_data['pvalue'].apply(lambda x: max(x, 1e-150))
        else:
            return "GWAS pvalue underflow."

    # Make sure all eQTLs are significant enough that this site is worth testing
    if min(eqtl_subset['pvalue']) > settings["eqtl_threshold"]:
        return "Insignificant eQTL top hit: -logp {0}".format(max([-math.log10(p) for p in eqtl_subset['pvalue']]))

    # Get MAFs from 1000 Genomes.
    # Filter out multi-allelic or non-polymorphic sites.
    # TODO: Make sure direction of MAFs is consistent between eQTL and GWAS
    # Currently, some MAFs may be > 0.5
    # (Could eventually be in separate function)
    # Get the region of interest from 1K genomes VCFs using tabix
    output = subprocess.check_output("tabix /mnt/lab_data/montgomery/shared/1KG/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz {0}:{1}-{2}".format(snp.chrom, snp.pos - window, snp.pos + window), shell=True).strip().split("\n")
    mafs = [[int(line.split('\t')[0]), int(line.split('\t')[1]), line.split('\t')[7]] for line in output if "MULTI_ALLELIC" not in line and ";AF=1;" not in line and ";AF=0;" not in line and "," not in line.split('\t')[4]]
    for m in mafs:
        m[2] = float(m[2].split(";")[1][3::])
    mafs = pd.DataFrame(mafs, columns=["chr_eqtl", "snp_pos", "Kgenomes_maf"])

   # TODO TODO: At this step filter out variants
    # whose same position appears twice or more. (can steal the code 
    # that is currently used in the FINEMAP pipeline only). This is important 
    # because otherwise at this level some sites with 2 SNPs might be merged into
    # 2x2 SNPs or something like that

    # TODO: Convert this to a join operation instead
    # Join the list of eQTL SNPs with the list of GWAS SNPs
    combined = pd.merge(gwas_data, eqtl_subset, on="snp_pos", suffixes=("_gwas", "_eqtl"))
    combined = pd.merge(combined, mafs, on=["snp_pos", "chr_eqtl"])

    # Filter out variants where 1K genomes MAF < 0.01. We can think more about
    # whether this is the best strategy, but for now it's best to do this, because
    # the overwhelming majority of variants with lower MAF end up getting filtered out
    # at later stages in the pipeline anyway.
    combined = combined[(combined['Kgenomes_maf'] > 0.01) & (combined['Kgenomes_maf'] < 0.99)]
 
    # For now, remove all positions that appear multiple times in the GWAS table.
    # This will avoid problems later in the pipeline, and doesn't remove too many SNPs anyway.
    dup_counts = {}
    for pos in combined['snp_pos']:
            dup_counts[pos] = dup_counts.get(pos, 0) + 1

    combined['dup_counts'] = [dup_counts[pos] for pos in combined['snp_pos']]
    combined = combined[combined['dup_counts'] == 1]

    # Check to make sure there are SNPs remaining; if not, just move on
    # to next gene.
    if combined.shape[0] == 0: 
        return "No overlapping SNPs in eQTL and GWAS"

    if min(combined['pvalue_gwas']) > settings["gwas_threshold"]:
        return "No significant GWAS SNPs are in eQTL dataset (too rare)"

    return combined

