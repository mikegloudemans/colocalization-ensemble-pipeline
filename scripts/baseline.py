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

    print locus.chrom
    print locus.pos

    combined = locus.data.copy()

    pval = 2
    pval2 = 1 # Slightly different baseline method; we only consider eQTLs directly overlapping lead GWAS SNP
    hit = combined[(combined["chr_gwas"] == locus.chrom) & (combined["snp_pos"] == locus.pos)]
    print hit.head()
    if hit.shape[0] == 1:
        # If the locus is in the data, return its p-value
        print hit["pvalue_eqtl"]
        print hit["pvalue_eqtl"].iloc[0]
        pval = hit["pvalue_eqtl"].iloc[0]
        pval2 = hit["pvalue_eqtl"].iloc[0]

    else:
        # Otherwise return best p-value with 10000 bp
        region = combined[(combined["chr_gwas"] == locus.chrom) & (abs(combined["snp_pos"] - locus.pos) <= 10000)]
        pval = min(region["pvalue_eqtl"])

    with open("{0}/{1}_baseline_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, locus.data.shape[0], locus.gwas_suffix, -1*math.log10(pval), -1*math.log10(pval2)))

    # Another option: If not in the data, return p-value of nearest LD buddy of SNPs that are in the data?

    # (Might even consider a few variations on this baseline model, like allowing another SNP or
    # allowing the best in the region...whichever one gives best results.)

    return pval
