# Enloc needs a separate dispatch script, because it 
# runs the whole genome at once

# So we need to run enloc and then we can pull out the results from
# the relevant files later from the main dispatch script

# -----------------------------------------------------------------

# Input:
#   - GWAS sumstats file (single trait is ideal, otherwise supply trait and filter down)
#   - eQTL sumstats file
#   - list of SNPs to test
#   - window size to consider

run_eqtl_dap = False

# Format eQTL files if needed, and run DAP
# Shouldn't be necessary every time
if run_eqtl_dap == True:
    subprocess.call("python enloc/run_dap_for_sims.py")
    pass

#
# Format GWAS file to include all loci of interest here...
# We'll just have to assume one trait/SNP pairing
#
# Would be nice to talk with William about this at some point to
# hear what his exact recommendation is.
#


