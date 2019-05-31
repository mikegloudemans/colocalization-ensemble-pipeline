# Author: Mike Gloudemans
# Author: Abhiram Rao

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
import os
import numpy as np

# Integration of code written by Ram
def run_gsmr(locus, window=1000000):
    
    subprocess.call("mkdir -p {4}/gsmr/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.tmpdir), shell=True)

    pval = 1

    dat = locus.data.copy()

    # We need beta and se for GWAS, or else this analysis won't work
    if "beta_gwas" not in dat or "se_gwas" not in dat:
        return "Need beta and SE for GWAS in order to run SMR."

    # Get rsids
    # TODO later: Move this to an earlier point so it doesn't have to be repeated
    # for plotting or anything
    
    # If both GWAS and eQTL have rsids, they'll have been renamed
    if "rsid_eqtl" in list(dat.columns.values):
        dat['rsid'] = dat['rsid_eqtl']

    if "rsid_gwas" in list(dat.columns.values):
        dat['rsid'] = dat['rsid_gwas']

    if "rsid" not in list(dat.columns.values):
        if "rsid_index_file" in locus.settings:

            # First, extract nearby variants using tabix
            stream = StringIO(subprocess.check_output("tabix {0} {1}:{2}-{3}".format(locus.settings["rsid_index_file"], locus.data['chr_gwas'][0], np.min(np.array(locus.data["snp_pos"])), np.max(np.array(locus.data["snp_pos"]))), shell=True))

            # For readability, load the header too
            # Load with pandas
            dbsnp = pd.read_csv(stream, sep="\t", header=None).iloc[:,:3]
            dbsnp = dbsnp.rename({0: "chr_gwas", 1: "snp_pos", 2:"rsid"}, axis="columns")
           
            dat = pd.merge(dbsnp, dat, left_on=["chr_gwas", "snp_pos"], right_on=["chr_gwas", "snp_pos"])

        else:
            return "No rsids specified; skipping locus."

    # See if we need to infer effect_af from the reference VCF
    if "effect_af_gwas" not in dat.columns.values or "effect_af_eqtl" not in dat.columns.values:
        # Code for this is currently in the COLOC script
        import coloc
        dat = coloc.get_mafs(locus, dat, window)
        if "effect_af_eqtl" not in dat:
            dat["effect_af_eqtl"] = dat["ref_af"]
        if "effect_af_gwas" not in dat:
            dat["effect_af_gwas"] = dat["ref_af"]

    dat.to_csv("{0}/gsmr/{1}/{2}_{3}/{4}/{5}_level{6}.csv".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), sep=",", index=False)

    # TODO: Display more precise error messages indicating why failed
    # tests have failed (probably due to not enough significant SNPs)
    try:
        pval = float(subprocess.check_output("Rscript run_gsmr.R {0}/gsmr/{1}/{2}_{3}/{4}/{5}_level{6} 2> /dev/null".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True))
    except:
        pass

    # Write results to the output file
    with open("{0}/{1}_gsmr_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, locus.data.shape[0], locus.gwas_suffix, -1*math.log10(pval)))

    return pval


