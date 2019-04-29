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
import gzip

import sys
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO


# Input: gwas file location, threshold of significance, minimum distance
# between two selected SNPs.  Output: A list of significant SNPs to test.
def select_test_snps_by_gwas(gwas_file, gwas_threshold, trait, settings, window=1000000):

    print("Selecting GWAS hits from {0}".format(gwas_file))

    with gzip.open(gwas_file) as f:
        if "debug" in settings and settings["debug"] == "True":
            gwas_table = pd.read_csv(f, sep="\t", nrows=500000, dtype=str)
        else:
            gwas_table = pd.read_csv(f, sep="\t", dtype=str)

    if trait == gwas_file.split("/")[-1]:
        subset = gwas_table[['chr', 'snp_pos', 'pvalue']].copy()
    else:
        subset = gwas_table[['chr', 'snp_pos', 'pvalue', 'trait']].copy()

    subset.loc[:,'pvalue'] = subset.loc[:,'pvalue'].astype(float)
    subset = subset[subset['pvalue'] <= gwas_threshold]
    if trait != gwas_file.split("/")[-1]:
        subset = subset[subset['trait'] == trait]

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


    # -1 here means unrestricted to a certain gene
    # TODO: Just make this part of the SNP object so that we can test SNP/gene pairs for equality
    snps_to_test = [(s, -1) for s in snps_to_test]

    return snps_to_test

# Input: gwas file location, threshold of significance, minimum distance
# between two selected SNPs.  Output: A list of significant SNPs to test.
def select_snps_from_list(list_file):

    gwas_table = pd.read_csv(list_file, sep="\s+", header=None)

    snp_list = zip(list(gwas_table.iloc[:,0]), list(gwas_table.iloc[:,1]), [-1]*gwas_table.shape[0])
    # Filter out X and Y for now
    snp_list = [s for s in snp_list if ("Y" not in s[0] and "X" not in s[0])]
    snp_list = [SNP.SNP(s) for s in snp_list]
    if gwas_table.shape[1] > 2:
        return zip(snp_list, list(gwas_table.iloc[:,2]))
    else:
        return zip(snp_list, [-1]*len(snp_list))



def select_test_snps_by_eqtl(eqtl_file, settings, subset_file=-1):

    print("Selecting eQTL hits from {0}".format(eqtl_file))

    snps_to_test = []

    # Load eQTL file
    with gzip.open(eqtl_file) as stream:

        # TODO: I may have tried this before, but I don't see what we don't
        # just convert chisq to p-value straight off the bat to make things simple?
        # But I shouldn't fix this until I'm about to run it on chisq-formatted samples again
        header = stream.readline().strip().split()
        if settings['eqtl_experiments'][eqtl_file]['eqtl_format'] == "chisq":
            assert 'chisq' in header
            mode = "chisq"
            pval_index = header.index("chisq")
            chisq_threshold = stats.chi2.isf(settings['selection_thresholds']['eqtl'],1)
        elif settings['eqtl_experiments'][eqtl_file]['eqtl_format'] == "effect_size":
            assert 'se' in header
            assert 'beta' in header
            se_index = header.index('se')
            beta_index = header.index('beta')
            mode = "effect_size"
        else:
            mode = "pvalue"
            pval_index = header.index("pvalue")

        if "feature" in header:
            feature_index = header.index("feature")
        else:
            feature_index = header.index("gene")
        pos_index = header.index("snp_pos")
        chrom_index = header.index("chr")

        gene_bests = {}

        if subset_file != -1:
            feature_subset = []
            with open(subset_file) as f:
                for line in f:
                    feature_subset.append(line.strip())
            feature_subset = set(feature_subset)

        i = 0
        for line in stream:
            i = i + 1

            # If debugging, only run a shorter chunk
            if i % 1000000 == 0:
                if "debug" in settings and settings["debug"] == True:
                    break
            data = line.split()
            feature = data[feature_index]
            if subset_file != -1:
                if feature not in feature_subset:
                    continue

            chrom = data[chrom_index].replace("chr", "")

            # Exclude sex chromosomes for now
            try:
                chrom = int(chrom)
            except:
                continue


            if mode == "chisq":
                chisq = float(data[pval_index])
                if chisq < chisq_threshold:
                    continue

                if feature not in gene_bests or gene_bests[feature][2] < chisq:
                    gene_bests[feature] = (chrom, data[pos_index], chisq)
            elif mode == "effect_size":
                try:
                    zscore = abs(float(data[beta_index]) / float(data[se_index]))
                    pvalue = stats.norm.sf(zscore) * 2
                except ValueError: # NA values
                    continue
                if feature not in gene_bests or gene_bests[feature][2] > pvalue:
                    gene_bests[feature] = (chrom, data[pos_index], pvalue)
            else:
                pvalue = float(data[pval_index])
                if pvalue > settings['selection_thresholds']['eqtl']:
                    continue

                if feature not in gene_bests or gene_bests[feature][2] > pvalue:
                    gene_bests[feature] = (chrom, data[pos_index], pvalue)

    # For each gene, determine the most significant SNP and add it to our list
    for gene in sorted(gene_bests.keys()):
        best = gene_bests[gene]
        if mode == "chisq":
            snps_to_test.append((SNP.SNP((best[0], best[1], stats.chi2.sf(best[2], 1))), gene))
        else:
            snps_to_test.append((SNP.SNP((best[0], best[1], best[2])), gene))

    return snps_to_test

# Load summary statistics for GWAS
def get_gwas_data(gwas_file, snp, settings, trait):

    window = settings["window"]

    # Get GWAS data using tabix
    #header = subprocess.check_output("zcat {0} | head -n 1".format(gwas_file), shell=True)

    with gzip.open(gwas_file, 'rb') as gwas:
        header = gwas.readline()

    raw_gwas = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(gwas_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True) + \
            subprocess.check_output("tabix {0} chr{1}:{2}-{3}".format(gwas_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True)
    gwas_table = pd.read_csv(StringIO(header + raw_gwas), sep="\t")
    gwas_table['pvalue'] = gwas_table['pvalue'].astype(float)

    if trait != gwas_file.split("/")[-1]:
        gwas_table = gwas_table[gwas_table["trait"] == trait]

    if gwas_table.shape[0] == 0:
        return "No GWAS summary statistics found at this locus."

    gwas_table['snp_pos'] = gwas_table['snp_pos'].astype(int)

    if "ref_allele_header" in settings['gwas_experiments'][gwas_file]:
        gwas_table['ref'] = gwas_table[settings['gwas_experiments'][gwas_file]['ref_allele_header']]
    if "alt_allele_header" in settings['gwas_experiments'][gwas_file]:
        gwas_table['alt'] = gwas_table[settings['gwas_experiments'][gwas_file]['alt_allele_header']]

    if "effect_allele" in list(gwas_table.columns.values):
        gwas_table['alt'] = gwas_table["effect_allele"]
        gwas_table['ref'] = gwas_table["non_effect_allele"]


    gwas_table['ref'] = gwas_table['ref'].apply(lambda x: x.upper())
    gwas_table['alt'] = gwas_table['alt'].apply(lambda x: x.upper())

    #
    # 'gwas_format' must be specified, to make sure the users know what they're doing.
    #
    # Possible settings for 'gwas_format':
    #   - case_control
    #   - effect_size
    #   - pval_only (requires effect direction)
    #

    if settings['gwas_experiments'][gwas_file]['gwas_format'] == 'case_control':
        assert 'log_or' in gwas_table and 'se' in gwas_table
        gwas_table['ZSCORE'] = gwas_table['log_or'] / gwas_table['se']
    elif settings['gwas_experiments'][gwas_file]['gwas_format'] == 'effect_size':
        if "log_or" in gwas_table:
            gwas_table['beta'] = gwas_table['log_or']
        assert 'beta' in gwas_table
        gwas_table['ZSCORE'] = gwas_table['beta'] / gwas_table['se']
        gwas_table['ZSCORE'] = gwas_table['ZSCORE'].fillna(0)
    elif settings['gwas_experiments'][gwas_file]['gwas_format'] == 'pval_only':
        assert 'pvalue' in gwas_table and ("effect_direction" in gwas_table or "direction" in gwas_table or "beta" in gwas_table)
        # Need to cap it at z-score of 40 for outrageous p-values (like with AMD / RPE stuff)
        if "effect_direction" in gwas_table:
            gwas_table['ZSCORE'] = pd.Series([min(x, 40) for x in stats.norm.isf(gwas_table["pvalue"] / 2)], index=gwas_table.index) * (2*(gwas_table["effect_direction"] == "+")-1)
        elif "beta" in gwas_table:
            # replace beta == NaN with 0 (beta == NaN for large p-values)
            gwas_table = gwas_table.fillna({'beta': 0})
            gwas_table['ZSCORE'] = pd.Series([min(x, 40) for x in stats.norm.isf(gwas_table["pvalue"] / 2)], index=gwas_table.index) * (2*(gwas_table["beta"] > 0)-1)
        else:
            gwas_table['ZSCORE'] = pd.Series([min(x, 40) for x in stats.norm.isf(gwas_table["pvalue"] / 2)], index=gwas_table.index) * (2*(gwas_table["direction"] == "+")-1)
    else:
        return "Improper GWAS format specification"

    return gwas_table

# Load summary statistics for eQTL
def get_eqtl_data(eqtl_file, snp, settings):

    window = settings["window"]

    # Get eQTL data using tabix
    with gzip.open(eqtl_file, 'rb') as eqtl:
        header = eqtl.readline()

    #header = subprocess.check_output("zcat {0} | head -n 1".format(eqtl_file), shell=True)

    raw_eqtls = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(eqtl_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True)
    raw_eqtls += subprocess.check_output("tabix {0} chr{1}:{2}-{3}".format(eqtl_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True)
    eqtls = pd.read_csv(StringIO(header + raw_eqtls), sep="\t", index_col=False)

    if eqtls.shape[0] == 0:
        return "Gene desert."

    # Should probably eventually remove the following four lines, just to simplify things
    '''
    if "ref_allele_header" in settings['eqtl_experiments'][eqtl_file]:
        eqtl_table['ref'] = eqtl_table[settings['eqtl_experiments'][eqtl_file]['ref_allele_header']]
    if "alt_allele_header" in settings['eqtl_experiments'][eqtl_file]:
        eqtl_table['alt'] = eqtl_table[settings['eqtl_experiments'][eqtl_file]['alt_allele_header']]
    '''

    eqtls['ref'] = eqtls['ref'].apply(lambda x: x.upper())
    eqtls['alt'] = eqtls['alt'].apply(lambda x: x.upper())

    eqtls['snp_pos'] = eqtls['snp_pos'].astype(int)#

    #
    # 'eqtl_format' must be specified, to make sure the users know what they're doing.
    #
    # Possible settings for 'eqtl_format':
    #   - tstat
    #   - effect_size
    #   - chisq
    #

    if settings['eqtl_experiments'][eqtl_file]['eqtl_format'] == 'tstat':
        assert 't-stat' or "tstat" in eqtls
        if "tstat" in eqtls:
            eqtls['t-stat'] = eqtls["tstat"]
        eqtls['ZSCORE'] = eqtls['t-stat']
    elif settings['eqtl_experiments'][eqtl_file]['eqtl_format'] == 'effect_size':
        assert 'beta' in eqtls
        eqtls['ZSCORE'] = eqtls['beta'] / eqtls['se']
        eqtls['pvalue'] = stats.norm.sf(abs(eqtls['beta'] / eqtls['se']))*2
    elif settings['eqtl_experiments'][eqtl_file]['eqtl_format'] == 'chisq':
        assert "chisq" in eqtls
        # Here we're dealing with RASQUAL data
        # where effect size is given by allelic imbalance percentage pi.
        # Use max function to protect against underflow in chi2 computation
        eqtls['pvalue'] = stats.chi2.sf(eqtls["chisq"],1)
        eqtls['ZSCORE'] = stats.norm.isf(eqtls['pvalue']/2) * (2 * (eqtls["pi"] > 0.5) - 1)
    elif settings['eqtl_experiments'][eqtl_file]['eqtl_format'] == 'pval_only':
        assert 'pvalue' in eqtls and ("effect_direction" in eqtls or "beta" in eqtls)
        # Need to cap it at z-score of 40 for outrageous p-values (like with AMD / RPE stuff)
        if "effect_direction" in eqtls:
            eqtls['ZSCORE'] = pd.Series([min(x, 40) for x in stats.norm.isf(eqtls["pvalue"] / 2)], index=eqtls.index) * (2*(eqtls["effect_direction"] == "+")-1)
        else:
            eqtls['ZSCORE'] = pd.Series([min(x, 40) for x in stats.norm.isf(eqtls["pvalue"] / 2)], index=eqtls.index) * (2*(eqtls["beta"] > 0)-1)

    else:
        return "Improper eQTL format specification"

    return eqtls

# Input: GWAS pandas dataframe, eQTL pandas dataframe, gene name as a string,
#   GWAS SNP as a tuple.
# Returns: a combined table of summary statistics, or None if we need to skip
#   the site due to insufficient data.
def combine_summary_statistics(gwas_data, eqtl_data, gene, snp, settings, unsafe=False, allow_insignificant_gwas=False):
    # TODO TODO TODO: Fix the SettingWithCopyWarning. It seems likely to be error-prone, according
    # to the manual!

    window = settings["window"]

    # Filter SNPs down to the gene of interest.
    eqtl_subset = eqtl_data[eqtl_data['gene'] == gene]

    if eqtl_subset.shape[0] == 0:
        return "The pre-specified gene was not tested in this sample."

    # Sometimes the GWAS SNP is outside of the range of eQTLs tested for a certain
    # gene, or on the outside fringe of the range. If this is the case, then skip it.
    # NOTE: Modify the 50000 cutoff if it doesn't seem like it's giving enough room for LD decay to fall off.

    if snp.pos > max(eqtl_subset['snp_pos']) - 50000 or snp.pos < min(eqtl_subset['snp_pos']) + 50000:
            return "SNP outside range."

    # If not explicitly allowing them, remove pvalues with danger
    # of underflow.

    # TODO: This should also be fixable

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
    if "screening_thresholds" in settings and "eqtl" in settings["screening_thresholds"]:
        if min(eqtl_subset['pvalue']) > settings["screening_thresholds"]["eqtl"]:
            return "Insignificant eQTL top hit: -logp {0}".format(max([-math.log10(p) for p in eqtl_subset['pvalue']]))

    # NOTE: At some point, we may want to do all 1K genomes filtering here by default, if necessary.

    combined = pd.merge(gwas_data, eqtl_subset, on="snp_pos", suffixes=("_gwas", "_eqtl"))
   
    # For now, remove all positions that appear multiple times in the GWAS table.
    # This will avoid problems later in the pipeline, and doesn't remove too many SNPs anyway.
    # NOTE: This might not actually be necessary at this point in the process, but we'll keep it just in case.
    dup_counts = {}
    for pos in combined['snp_pos']:
            dup_counts[pos] = dup_counts.get(pos, 0) + 1

    combined['dup_counts'] = [dup_counts[pos] for pos in combined['snp_pos']]
    combined = combined[combined['dup_counts'] == 1]

    # Check to make sure there are SNPs remaining; if not, just move on
    # to next gene.
    if combined.shape[0] == 0:
        return "No overlapping SNPs in eQTL and GWAS"

    # Check to make sure we still have significant GWAS hits and eQTLs, if desired
    if "screening_thresholds" in settings and "gwas" in settings["screening_thresholds"]:
        if min(combined['pvalue_gwas']) > settings["screening_thresholds"]["gwas"]:
            return "No significant GWAS SNPs are in eQTL dataset (too rare)"

    return combined

