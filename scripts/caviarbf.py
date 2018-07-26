# Author: Mike Gloudemans
#
# Run eCAVIAR algorithm for colocalization.
# Store full results in the output directory for this run.
# Return the CLPP score.
#

import subprocess
from scipy import stats
import finemap
import sys

def run_caviarbf(locus, window=500000):

    # Run the part shared between this and the finemap
    # pipeline, but only if the finemap pipeline hasn't already
    # run this part.
    if "finemap" not in locus.settings["methods"] and "ecaviar" not in locus.settings["methods"]:
        finemap.prep_finemap(locus, window)

    # Run CAVIARBF
    command = '/users/mgloud/software/caviarbf/caviarbf/caviarbf \
              -o {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_caviarbf_results_gwas \
              -r {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.fixed.ld \
              -z {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z \
              -t 0 \
              -a 1 \
              --appr \
              -n {7} \
              -c 1'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir, locus.settings["gwas_experiments"][locus.gwas_file]["N"])
    subprocess.check_call(command, shell=True)

    command = '/users/mgloud/software/caviarbf/caviarbf/caviarbf \
              -o {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_caviarbf_results_eqtl \
              -r {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.fixed.ld \
              -z {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z \
              -t 0 \
              -a 1 \
              --appr \
              -n {7} \
              -c 1'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir, locus.settings["eqtl_experiments"][locus.eqtl_file]["N"])
    subprocess.check_call(command, shell=True)

    sys.exit()

    # Parse eCAVIAR results to compute CLPP score
    command = "awk '{{sum += $2}} END {{print sum}}' ../tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_ecaviar_results_col".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level) 
    clpp = float(subprocess.check_output(command, shell=True))
    
    snps_tested = int(check_output("wc -l ../tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_ecaviar_results_col".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True))

    # Add results to the desired file
    with open("{0}/{1}_ecaviar_clpp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, snps_tested, clpp))

    return clpp
