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

def run_rtc(locus, window=500000):

    # Directions for running RTC are online at https://qtltools.github.io/qtltools/pages/mode_rtc.html

    # Step 1: Run the permutation pass

    subprocess.call("mkdir -p {0}/rtc".format(locus.tmpdir), shell=True)

    # (hopefully it will still run without any covariates given?)
    subprocess.check_call("QTLtools cis --vcf {0} --bed {1} --permute 200 --out {2}/rtc/permutations_all.txt".format(locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_genos"], locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_phenos"], locus.tmpdir), shell=True)
    subprocess.check_call("gzip -c {0}/rtc/permutations_all.txt".format(locus.tmpdir), shell=True)

    # Run FDR correction
    # I'm not sure where this script comes from? Maybe is in QTLtools directory?
    subprocess.check_call("Rscript ./script/runFDR_cis.R {0}/rtc/permutations_all.txt.gz 0.05 {0}/rtc/permutations_all".format(locus.tmpdir))

    # Step 2: Run RTC

    # Need to create GWAS file
    with open("{0}/rtc/gwas.txt".format(locus.tmpdir), "w") as w:
        w.write("{0}\t{1}\n".format(locus[0], locus[1]))
    # That was easy, I hope

    # NOTE: This step can be tuned with a "--conditional" option to run for independent conditional eQTLs (see website for more info)
    subprocess.check_call("QTLtools rtc --vcf {0} --bed {1} --cov --hotspot /users/mgloud/projects/brain_gwas/data/rtc/hotspots_b37_hg19.bed --gwas-cis {2}/rtc/gwas.txt {2}/rtc/permutations_all.significant.txt --normal --out {2}/rtc/rtc_results.txt".format(locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_genos"], locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_phenos"], locus.tmpdir), shell=True)
    
    # Parse RTC results to compute CLPP score
    # TODO

    # Add results to the desired file
    #with open("{0}/{1}_ecaviar_clpp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
    #    a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, snps_tested, clpp))

    sys.exit()
    return rtc
