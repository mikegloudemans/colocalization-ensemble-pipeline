# Author: Mike Gloudemans
#
# Run Giambortolomei et al. algorithm for colocalization.
# Store full results in the output directory for this run.
# Return a tuple of (H0, H1, H2, H3, H4) posterior probabilities.
#

import subprocess
from scipy import stats
from shutil import copyfile

def run_coloc(locus, window=500000):

    prep_coloc(locus, window)
    launch_coloc(locus, window)

def prep_coloc(locus, window):
    # Write a simple table-formatted file containing the necessary information
    # for coloc to run

    conditional_level = 0

    # TODO: Standardize the input format so this can be used on all input files.
    # TODO: Test this on brain data instead, so we can have a more appropriate input.
        # The RPE stuff isn't really right for this coloc program.

    subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)

    locus.data.to_csv("/users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, conditional_level), sep="\t", columns=["pvalue_gwas", "pvalue_eqtl", "allele_freq", "rsid", "Ncases", "Ncontrols"], index=False)

def launch_coloc(locus, window):
    # Launch an R script to run the coloc package.

    # Parse results

    # PURGE tmp files after that - as simple as removing the whole SNP directory.
    pass

