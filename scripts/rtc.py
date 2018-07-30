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
import os

def run_rtc(locus, window=500000):

    # Directions for running RTC are online at https://qtltools.github.io/qtltools/pages/mode_rtc.html

    # Step 1: Run the permutation pass

    subprocess.call("mkdir -p {0}/rtc".format(locus.tmpdir), shell=True)
    # This next line is important; otherwise if one run falls through,
    # its results will be counted for the next run
    subprocess.call("rm {0}/rtc/rtc_results.txt 2> /dev/null".format(locus.tmpdir), shell=True)

    # (hopefully it will still run without any covariates given?)
    print "QTLtools cis --vcf {0} --bed {1} --permute 200 --out {2}/rtc/permutations_all.txt".format(locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_genos"], locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_phenos"], locus.tmpdir)
    subprocess.check_call("QTLtools cis --vcf {0} --bed {1} --permute 200 --out {2}/rtc/permutations_all.txt".format(locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_genos"], locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_phenos"], locus.tmpdir), shell=True)

    # TODO: Implement this FDR filtering once we get the simulations going for
    # genome-wide simulations. Then we can consider the exact level of FDR filtering
    # to be an experimental parameter.
    # Run FDR correction
    # NOTE: RTC is really made to run across the whole genome at once after calling eQTLs...
    # consider modifying the pipeline eventually so that it can use a whole "fake genome" to
    # calculate FDR, in fairness to the method.
    #print "Rscript /users/mgloud/software/qtltools/script/runFDR_cis.R {0}/rtc/permutations_all.txt.gz 0.05 {0}/rtc/permutations_all".format(locus.tmpdir)
    #subprocess.check_call("Rscript /users/mgloud/software/qtltools/script/runFDR_cis.R {0}/rtc/permutations_all.txt.gz 0.05 {0}/rtc/permutations_all".format(locus.tmpdir))


    # Step 2: Run RTC

    # Need to create GWAS file
    # TODO: Make sure we're feeding RTC the top simulated GWAS SNP instead of the 
    # one that is known by the simulator but isn't known a priori to an observer
    with open("{0}/rtc/gwas.txt".format(locus.tmpdir), "w") as w:
        w.write("{0}_{1}\n".format(locus.chrom, locus.pos))
    # That was easy, I hope

    # NOTE: This step can be tuned with a "--conditional" option to run for independent conditional eQTLs (see website for more info)
    print "QTLtools rtc --vcf {0} --bed {1} --hotspots /users/mgloud/projects/brain_gwas/data/rtc/hotspots_b37_hg19.bed --gwas-cis {2}/rtc/gwas.txt {2}/rtc/permutations_all.txt --normal --out {2}/rtc/rtc_results.txt".format(locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_genos"], locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_phenos"], locus.tmpdir)
    subprocess.check_call("QTLtools rtc --vcf {0} --bed {1} --hotspots /users/mgloud/projects/brain_gwas/data/rtc/hotspots_b37_hg19.bed --gwas-cis {2}/rtc/gwas.txt {2}/rtc/permutations_all.txt --normal --out {2}/rtc/rtc_results.txt".format(locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_genos"], locus.settings["eqtl_experiments"][locus.eqtl_file]["rtc_phenos"], locus.tmpdir), shell=True)
   
    # NOTE: A lot of the time, RTC won't run at all if the top eQTL SNP and the top
    # GWAS SNP are separated by one of the recombination hotspots. This is probably
    # okay; it seems like a feature of this particular method if anything.

    # Parse RTC results to compute CLPP score
    # TODO
    with open("{0}/rtc/rtc_results.txt".format(locus.tmpdir)) as f:
        f.readline()
        rtc = f.readline()
        if rtc=="":
            rtc = 0
        else:
            rtc = float(rtc.strip().split()[19])

    # Add results to the desired file
    with open("{0}/{1}_rtc_score_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
                a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, rtc, locus.gwas_suffix))
    
    return rtc
