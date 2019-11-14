import glob
import subprocess
import gzip

# Get sim eQTL files
eqtl_files = glob.glob("/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/eqtl/eqtl_sumstats*.gz")

# Read in eQTL files and write them all to a separate file
for i in range(len(eqtl_files)):
    print i
    with open("/users/mgloud/projects/brain_gwas/scripts/enloc/tmp/eqtl_pre_dap.txt", "w") as w:
        w.write("SNP\tlocus\tz-val\n")
        with gzip.open(eqtl_files[i]) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                feature = data[0]
                chrom = data[4]
                pos = data[5]
                zscore = data[11]
                w.write("{0}:{1}\tLoc{2}\t{3}\n".format(chrom, pos, i, zscore))

    # Run DAP-G for eQTLs
    subprocess.call("/users/mgloud/software/enloc/integrative/src/dap1/dap1 -d /users/mgloud/projects/brain_gwas/scripts/enloc/tmp/eqtl_pre_dap.txt > /users/mgloud/projects/brain_gwas/scripts/enloc/tmp/sim_eqtls/eqtl{0}_dap.txt".format(i), shell=True)


