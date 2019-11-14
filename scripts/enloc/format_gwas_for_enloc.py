import glob
import gzip
import subprocess

# Input should be a set of GWAS files with the desired data...
gwas = glob.glob("/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/gwas/*.gz")

subprocess.check_call("rm -f /users/mgloud/projects/brain_gwas/scripts/enloc/tmp/gwas_data.gz", shell=True)

# Read in GWAS files and write them all to a separate file
with gzip.open("/users/mgloud/projects/brain_gwas/scripts/enloc/tmp/gwas_data.gz", "a") as a:
    for i in range(len(gwas)):
        with gzip.open(gwas[i]) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                chrom = data[3]
                pos = data[4]
                zscore = data[10]
                a.write("chr{0}:{1}\tLoc{2}\t{3}\n".format(chrom, pos, i, zscore))

# Run perl utilities to generate required info to run enloc
subprocess.call("perl /users/mgloud/software/enloc/integrative/dev/utility/get_max_pip.pl /users/mgloud/projects/brain_gwas/scripts/enloc/tmp/sim_eqtls > /users/mgloud/projects/brain_gwas/scripts/enloc/tmp/pip.info", shell=True)
subprocess.call("perl /users/mgloud/software/enloc/integrative/dev/utility/summarize_dapg_sig.pl /users/mgloud/projects/brain_gwas/scripts/enloc/tmp/sim_eqtls > /users/mgloud/projects/brain_gwas/scripts/enloc/tmp/eqtl.sig.summary", shell=True) 

# Write parameter file for enloc
