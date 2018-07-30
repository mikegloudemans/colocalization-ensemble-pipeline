# Author: Mike Gloudemans
#
# Run eCAVIAR algorithm for colocalization.
# Store full results in the output directory for this run.
# Return the CLPP score.
#

import subprocess
from scipy import stats
from scipy import misc
import finemap
import sys
import math

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
              -n {7} \
              -c 1'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir, locus.settings["gwas_experiments"][locus.gwas_file]["N"])
    print command
    subprocess.check_call(command, shell=True)


    command = '/users/mgloud/software/caviarbf/caviarbf/caviarbf \
              -o {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_caviarbf_results_eqtl \
              -r {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.fixed.ld \
              -z {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z \
              -t 0 \
              -a 1 \
              -n {7} \
              -c 1'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir, locus.settings["eqtl_experiments"][locus.eqtl_file]["N"])
    subprocess.check_call(command, shell=True)

    # Parse eCAVIAR results to compute CLPP score
    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_caviarbf_results_gwas".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)) as f:
        f.readline()
        f.readline()
        bf = []
        for line in f:
            bf.append(float(line.strip().split()[0]))
    denom = misc.logsumexp(bf)
    gwas_bf = [b - denom for b in bf]

    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_caviarbf_results_eqtl".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)) as f:
        f.readline()
        f.readline()
        bf = []
        for line in f:
            bf.append(float(line.strip().split()[0]))
    denom = misc.logsumexp(bf)
    eqtl_bf = [b - denom for b in bf]

    assert len(eqtl_bf) == len(gwas_bf)

    coloc_probs = [math.exp(eqtl_bf[i]) * math.exp(gwas_bf[i]) for i in range(len(gwas_bf))]
    clpp = sum(coloc_probs)

    snps_tested = int(subprocess.check_output("wc -l {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_caviarbf_results_eqtl".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), shell=True).split()[0]) - 2

    # Add results to the desired file
    with open("{0}/{1}_caviarbf_clpp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, snps_tested, clpp))

    print clpp
    return clpp
