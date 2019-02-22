# Author: Mike Gloudemans
#
# Run Giambortolomei et al. algorithm for colocalization.
# Store full results in the output directory for this run.
# Return a tuple of (H0, H1, H2, H3, H4) posterior probabilities.
#

# NOTE: Most of the time spent here currently is just from
# reloading the R packages for every single variant. Probably
# not worth changing the workflow though, since it still is at least
# as fast as FINEMAP and other algorithms.
#
# NOTE: That probably won't be true anymore once we start loading in MAFs
#

import subprocess
from scipy import stats
from shutil import copyfile
import sys
if sys.version_info[0] < 3: 
   from StringIO import StringIO
else:
   from io import StringIO
import pandas as pd

def run_coloc(locus, window=500000):

    prep_coloc(locus, window)
    return launch_coloc(locus, window)

def prep_coloc(locus, window):
    # Write a simple table-formatted file containing the necessary information
    # for coloc to run
    subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3} > /dev/null".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
    data = locus.data.copy()
    data = get_mafs(locus, data, window)

    if "effect_af_eqtl" not in data:
        data["effect_af_eqtl"] = data["ref_af"]
    if "effect_af_gwas" not in data:
        data["effect_af_gwas"] = data["ref_af"]

    data.to_csv("/users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), sep="\t", columns=["pvalue_gwas", "pvalue_eqtl", "effect_af_eqtl", "effect_af_gwas", "ZSCORE_gwas", "ZSCORE_eqtl"], index=False)

def launch_coloc(locus, window):

    # For quantitative traits there is no "s";
    # For case-control there may be no "beta" or "varbeta"
    if "s" not in locus.settings["gwas_experiments"][locus.gwas_file]:
        locus.settings["gwas_experiments"][locus.gwas_file]["s"] = 0

    # Launch an R script to run the coloc package.
    coloc_prob_h4 = float(subprocess.check_output("Rscript /users/mgloud/projects/brain_gwas/scripts/run_coloc.R /users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv {6} {7} {8} {9}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, \
        locus.settings["gwas_experiments"][locus.gwas_file]["N"], \
        locus.settings["gwas_experiments"][locus.gwas_file]["s"], \
        locus.settings["gwas_experiments"][locus.gwas_file]["type"], \
        locus.settings["eqtl_experiments"][locus.eqtl_file]["N"]), shell=True))
    # ^ Format of command line arguments: (N-gwas s-gwas type-gwas N-eqtl)

    num_sites = int(subprocess.check_output("wc -l /users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True).split()[0])-1

    # Write COLOC results to the desired file.
    with open("{0}/{1}_coloc_h4pp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, num_sites, coloc_prob_h4, locus.gwas_suffix))

    # Purge tmp files
    subprocess.call("rm -f /users/mgloud/projects/brain_gwas/tmp/coloc/{0}/{1}_{2}/{3}/{4}_level{5}.csv".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True)

    print coloc_prob_h4

    return coloc_prob_h4

# Get MAFs for data based on a reference VCF file
# Eliminate any data points for which MAFs can't be found
def get_mafs(locus, combined, window):
   
    ref = locus.settings["ref_genomes"][locus.settings["eqtl_experiments"][locus.eqtl_file]["ref"]]
    filename = ref["file"].format(locus.chrom)
    
    # First, extract nearby variants using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 500 | grep \\#CHROM | cut -f1,2,3,4,5,8".format(filename), shell=True).strip().split()
    if "chr_prefix" in ref and ref["chr_prefix"] == "chr":
        stream = StringIO(subprocess.check_output("tabix {8} {1}:{6}-{7} | cut -f1,2,3,4,5,8".format(locus.gwas_suffix, "chr" + str(locus.chrom), locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.pos-window, locus.pos+window, filename), shell=True))
    else:
        stream = StringIO(subprocess.check_output("tabix {8} {1}:{6}-{7} | cut -f1,2,3,4,5,8".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.pos-window, locus.pos+window, filename), shell=True))

    # For readability, load the header too
    # Load with pandas
    vcf = pd.read_csv(stream, sep=r"\s*", names=header)

    # Remove variants not in the GWAS table
    vcf["POS"] = (vcf["POS"]).astype(int)
    vcf = vcf[vcf["POS"].isin(list(combined["snp_pos"]))]

    # Remove variants with position appearing multiple times
    dup_counts = {}
    for pos in vcf["POS"]:
            dup_counts[pos] = dup_counts.get(pos, 0) + 1
    vcf["dup_counts"] = [dup_counts[pos] for pos in vcf['POS']]
    vcf = vcf[vcf["dup_counts"] == 1]

    # Remove multiallelic variants with only one entry in VCF
    l = lambda x: "," not in x
    vcf = vcf[vcf["REF"].apply(l) & vcf["ALT"].apply(l)]

    # Remove monoallelic variants.
    # Allele frequency might be input as counts or as percentages,
    # so handle this.
    if "af_attribute" in ref:
        af_id = ref["af_attribute"]
        def fn(x):
            info = [s for s in x.split(";") if s.startswith(af_id + "=")][0]
            af = float(info.split("=")[1])
            return af
        vcf['ref_af'] = vcf["INFO"].apply(fn)
        vcf = vcf[(vcf['ref_af'] > 0.01) & (1-vcf['ref_af'] > 0.01)]
    else:
        ac_id = ref["ac_attribute"]
        an = 2*ref["N"] # Assume 2 alleles per person
        def fn(x):
            info = [s for s in x.split(";") if s.startswith(ac_id + "=")][0]
            ac = float(info.split("=")[1])
            af = ac*1.0/an
            return af
        vcf['ref_af'] = vcf["INFO"].apply(fn)
        vcf = vcf[(vcf['ref_af'] > 0.01) & (1-vcf['ref_af'] > 0.01)]

    if "POS" in list(combined.columns.values):
        combined = combined.drop(columns=['POS'])

    # Remove variants where alt/ref don't match between GWAS/eQTL and VCF
    # Flipped is okay. A/C and C/A are fine, A/C and A/G not fine.
    merged = pd.merge(combined, vcf, left_on="snp_pos", right_on="POS")

    keep_indices = \
            (((merged['ref_gwas'] == merged['REF']) & (merged['alt_gwas'] == merged['ALT'])) | \
            ((merged['alt_gwas'] == merged['REF']) & (merged['ref_gwas'] == merged['ALT']))) & \
            (((merged['ref_eqtl'] == merged['REF']) & (merged['alt_eqtl'] == merged['ALT'])) | \
            ((merged['alt_eqtl'] == merged['REF']) & (merged['ref_eqtl'] == merged['ALT'])))

    keep = merged['POS'][keep_indices]
    merged = merged[merged['POS'].isin(list(keep))]

    return merged
