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
def select_test_snps_by_gwas(gwas_file, gwas_threshold, window=1000000):

    print("Selecting GWAS hits from {0}".format(gwas_file))

    stream = StringIO(subprocess.check_output("zcat {0}".format(gwas_file), shell=True))
    gwas_table = pd.read_csv(stream, sep="\t")
    subset = gwas_table[['chr', 'snp_pos', 'pvalue']]
    # TODO: Fix this line! Something is wrong with it I guess
    subset['pvalue'] = subset['pvalue'].astype(float)
    subset = subset[subset['pvalue'] <= gwas_threshold]

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
            chisq_threshold = stats.chi2.isf(settings['eqtl_threshold'],1)
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
            if i % 3000000 == 0:
                # NOTE temporary!
                pass
                #break
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
            else:
                pvalue = float(data[pval_index])
                if pvalue > settings['eqtl_threshold']:
                    continue

                if feature not in gene_bests or gene_bests[feature][2] > pvalue:
                    gene_bests[feature] = (chrom, data[pos_index], pvalue)

    # For each gene, determine the most significant SNP and add it to our list
    for gene in gene_bests:
        best = gene_bests[gene]
        if mode == "chisq":
            snps_to_test.append((SNP.SNP((best[0], best[1], stats.chi2.sf(best[2], 1))), gene))
        else:
            snps_to_test.append((SNP.SNP((best[0], best[1], best[2])), gene))

    return snps_to_test

# Load summary statistics for GWAS
def get_gwas_data(gwas_file, snp, settings):

    window = settings["window"]

    # Get GWAS data using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 1".format(gwas_file), shell=True)

    raw_gwas = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(gwas_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True) + \
            subprocess.check_output("tabix {0} chr{1}:{2}-{3}".format(gwas_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True)
    gwas_table = pd.read_csv(StringIO(header + raw_gwas), sep="\t")

    if gwas_table.shape[0] == 0:
        return "No GWAS summary statistics found as this locus."

    gwas_table['snp_pos'] = gwas_table['snp_pos'].astype(int)
    
    if "ref_allele_header" in settings['gwas_experiments'][gwas_file]:
        gwas_table['ref'] = gwas_table[settings['gwas_experiments'][gwas_file]['ref_allele_header']]
    if "alt_allele_header" in settings['gwas_experiments'][gwas_file]:
        gwas_table['alt'] = gwas_table[settings['gwas_experiments'][gwas_file]['alt_allele_header']]

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
        assert 'beta' in gwas_table
        gwas_table['ZSCORE'] = gwas_table['beta'] / gwas_table['se']
    elif settings['gwas_experiments'][gwas_file]['gwas_format'] == 'pval_only':
        assert 'pvalue' in gwas_table #and "direction" in gwas_table
        # Need to cap it at z-score of 40 for outrageous p-values (like with AMD / RPE stuff)
        gwas_table['ZSCORE'] = pd.Series([min(x, 40) for x in stats.norm.isf(gwas_table["pvalue"] / 2)]) # * (2*(gwas_table["direction"] == "+")-1)
    else:
        return "Improper GWAS format specification"

    return gwas_table

# Load summary statistics for eQTL
def get_eqtl_data(eqtl_file, snp, settings):

    window = settings["window"]

    # Get eQTL data using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 1".format(eqtl_file), shell=True)
    raw_eqtls = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(eqtl_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True)
    eqtls = pd.read_csv(StringIO(header + raw_eqtls), sep="\t")

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
    elif settings['eqtl_experiments'][eqtl_file]['eqtl_format'] == 'chisq':
        assert "chisq" in eqtls
        # Here we're dealing with RASQUAL data
        # where effect size is given by allelic imbalance percentage pi.
        # Use max function to protect against underflow in chi2 computation
        eqtls['pvalue'] = stats.chi2.sf(eqtls["chisq"],1)
        eqtls['ZSCORE'] = stats.norm.isf(eqtls['pvalue']/2) * (2 * (eqtls["pi"] > 0.5) - 1)
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
    # TODO: Figure out a better way of specifying how to do this thresholding; which ones to test
    # to make things fast
    #if min(eqtl_subset['pvalue']) > settings["eqtl_threshold"]:
    #    return "Insignificant eQTL top hit: -logp {0}".format(max([-math.log10(p) for p in eqtl_subset['pvalue']]))

    # Get MAFs from 1000 Genomes.
    # Filter out multi-allelic or non-polymorphic sites.
    # TODO TODO TODO: Make sure direction of MAFs is consistent between eQTL and GWAS
    # It only works right now because we're only using MAF for filtering
    # Currently, some MAFs may be > 0.5
    # (Could eventually be in separate function)
    # Get the region of interest from 1K genomes VCFs using tabix
    #output = subprocess.check_output("tabix /mnt/lab_data/montgomery/shared/1KG/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz {0}:{1}-{2}".format(snp.chrom, snp.pos - window, snp.pos + window), shell=True).strip().split("\n")
    #mafs = [[int(line.split('\t')[0]), int(line.split('\t')[1]), line.split('\t')[7]] for line in output if "MULTI_ALLELIC" not in line and ";AF=1;" not in line and ";AF=0;" not in line and "," not in line.split('\t')[4]]
    #for m in mafs:
    #    m[2] = float(m[2].split(";")[1][3::])
    #mafs = pd.DataFrame(mafs, columns=["chr_eqtl", "snp_pos", "Kgenomes_maf"])

    # TODO TODO: At this step filter out variants
    # whose same position appears twice or more. (can steal the code 
    # that is currently used in the FINEMAP pipeline only). This is important 
    # because otherwise at this level some sites with 2 SNPs might be merged into
    # 2x2 SNPs or something like that

    # TODO: Convert this to an indexed join operation instead
    # Join the list of eQTL SNPs with the list of GWAS SNPs
    combined = pd.merge(gwas_data, eqtl_subset, on="snp_pos", suffixes=("_gwas", "_eqtl"))
    #combined = pd.merge(combined, mafs, on=["snp_pos", "chr_eqtl"])

    # Filter out variants where 1K genomes MAF < 0.01. We can think more about
    # whether this is the best strategy, but for now it's best to do this, because
    # the overwhelming majority of variants with lower MAF end up getting filtered out
    # at later stages in the pipeline anyway.
    #combined = combined[(combined['Kgenomes_maf'] > 0.01) & (combined['Kgenomes_maf'] < 0.99)]
 
    # NOTE: I don't think we really need to do this anymore; however, we DO want to make sure
    # the same exact variant/rsid doesn't appear twice in the input. That would be a malformed
    # input, but it's happened before nevertheless.
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

    # TODO: Figure out how this works out when we're filtering based on eQTLs...
    # as it is right now, there's potential for trouble
    #if not allow_insignificant_gwas and min(combined['pvalue_gwas']) > settings["gwas_threshold"]:
    #    return "No significant GWAS SNPs are in eQTL dataset (too rare)"

    return combined

