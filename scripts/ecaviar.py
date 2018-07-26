# Author: Mike Gloudemans
#
# Run eCAVIAR algorithm for colocalization.
# Store full results in the output directory for this run.
# Return the CLPP score.
#

import subprocess
from scipy import stats
import finemap

def run_ecaviar(locus, window=500000):

    # Run the part shared between this and the finemap
    # pipeline, but only if the finemap pipeline hasn't already
    # run this part.
    if "finemap" not in locus.settings["methods"]:
        finemap.prep_finemap(locus, window)

    # Run eCAVIAR
    # NOTE: Redirecting output to /dev/null because eCAVIAR prints the entire
    # LD matrix, which is just bad when you have thousands of loci. Fix this if
    # you need to see eCAVIAR output.
    command = '/srv/persistent/bliu2/tools/caviar/CAVIAR-C++/eCAVIAR \
              -o {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_ecaviar_results \
              -l {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.fixed.ld \
              -l {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.fixed.ld \
              -z {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z \
              -z {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z \
              -c 1 \
              -r 0.95 > /dev/null'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)
    subprocess.check_call(command, shell=True)

    # Parse eCAVIAR results to compute CLPP score
    command = "awk '{{sum += $2}} END {{print sum}}' {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_ecaviar_results_col".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir) 
    clpp = float(subprocess.check_output(command, shell=True))

    snps_tested = int(check_output("wc -l {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_ecaviar_results_col".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), shell=True))

    # Add results to the desired file
    with open("{0}/{1}_ecaviar_clpp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, snps_tested, clpp))

    return clpp
