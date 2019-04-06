# Author: Mike Gloudemans

import subprocess
from scipy import stats
from shutil import copyfile
import sys
if sys.version_info[0] < 3: 
   from StringIO import StringIO
else:
   from io import StringIO
import pandas as pd
import math

def run_baseline(locus):
    
    combined = locus.data.copy()

    pval = 2
    pval2 = 1 # Slightly different baseline method; we only consider eQTLs directly overlapping lead GWAS SNP
    pval3 = 2 # Will combine p-values from the GWAS and eQTL
    hit = combined[(combined["chr_gwas"] == locus.chrom) & (combined["snp_pos"] == locus.pos)]
    print hit.head()
    if hit.shape[0] == 1:
        # If the locus is in the data, return its p-value
        pval = hit["pvalue_eqtl"].iloc[0]
        pval2 = hit["pvalue_eqtl"].iloc[0]
        pval3 = hit["pvalue_eqtl"].iloc[0] * hit["pvalue_gwas"].iloc[0]

    else:
        # Otherwise return best p-value within 10000 bp
        region = combined[(combined["chr_gwas"] == locus.chrom) & (abs(combined["snp_pos"] - locus.pos) <= 10000)]
        pval = min(region["pvalue_eqtl"])
        i = [i for enumerate(i, p) in region["pvalue_eqtl"] if (math.log(p) - math.log(pval < 0.0000001)][0]
        # For combined GWAS and p-value score, get the GWAS pvalue at that SNP instead of the root one
        gwas_p = region["pvalue_gwas"][i]

        pval3 = pval * gwas_p


    with open("{0}/{1}_baseline_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, locus.data.shape[0], locus.gwas_suffix, -1*math.log10(pval), -1*math.log10(pval2), -1*math.log10(pval3+1e-150)))

    # Another option: If not in the data, return p-value of nearest LD buddy of SNPs that are in the data?

    # (Might even consider a few variations on this baseline model, like allowing another SNP or
    # allowing the best in the region...whichever one gives best results.)

    return pval
