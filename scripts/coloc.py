# Author: Mike Gloudemans
#
# Run Giambortolomei et al. algorithm for colocalization.
# Store full results in the output directory for this run.
# Return a tuple of (H0, H1, H2, H3, H4) posterior probabilities.
#

# NOTE: Most of the time spent here currently is just from
# reloading the R packagaes for every single variant. Probably
# not worth changing the workflow though, since it still is at least
# as fast as FINEMAP and other algorithms.

import subprocess
from scipy import stats
from shutil import copyfile

def run_coloc(locus, window=500000):

    prep_coloc(locus, window)
    return launch_coloc(locus, window)

def prep_coloc(locus, window):
    # Write a simple table-formatted file containing the necessary information
    # for coloc to run

    # TODO: Standardize the input format so this can be used on all input files.

    subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)

    #locus.data.to_csv("/users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), sep="\t", columns=["pvalue_gwas", "pvalue_eqtl", "allele_freq"], index=False)
    locus.data.to_csv("/users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), sep="\t", columns=["pvalue_gwas", "pvalue_eqtl", "Kgenomes_maf"], index=False)

def launch_coloc(locus, window):
    # Launch an R script to run the coloc package.
    # NOTE: UNSAFE call here. Fix these all over the code. I should make a safe subprocess module that makes it impossible to accidentally remove the hard drive or whatever
    coloc_prob_h4 = float(subprocess.check_output("Rscript /users/mgloud/projects/brain_gwas/scripts/run_coloc.R /users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv {6} {7} {8} {9}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, \
        locus.settings["gwas_experiments"][locus.gwas_file]["N"], \
        locus.settings["gwas_experiments"][locus.gwas_file]["s"], \
        locus.settings["gwas_experiments"][locus.gwas_file]["type"], \
        locus.settings["eqtl_experiments"][locus.eqtl_file]["N"]), shell=True))
    # ^ Format of command line arguments: (N-gwas s-gwas type-gwas N-eqtl)

    num_sites = int(subprocess.check_output("wc -l /users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True).split()[0])-1

    # Write COLOC results to the desired file.
    with open("{0}/{1}_coloc_h4pp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
            a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, num_sites, coloc_prob_h4))

    # Purge tmp files
    subprocess.call("rm -f /users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True)

    return coloc_prob_h4
