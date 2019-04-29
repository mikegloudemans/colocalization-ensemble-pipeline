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
import gzip


def run_coloc(locus, window=500000):
    prep_coloc(locus, window)
    return launch_coloc(locus, window)

def prep_coloc(locus, window):
    # Write a simple table-formatted file containing the necessary information
    # for coloc to run
    subprocess.call("mkdir -p {0}/coloc/{1}/{2}_{3}/{4} > /dev/null".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
    data = locus.data.copy()
    if "effect_af_eqtl" not in data or "effect_af_gwas" not in data:
        data = get_mafs(locus, data, window)

        if "effect_af_eqtl" not in data:
            data["effect_af_eqtl"] = data["ref_af"]
        if "effect_af_gwas" not in data:
            data["effect_af_gwas"] = data["ref_af"]

    #data = data.dropna(axis=0) # remove incomplete rows to prevent error in run_coloc.R
    data.to_csv("{0}/coloc/{1}/{2}_{3}/{4}/{5}_level{6}.csv".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), sep=",", columns=["pvalue_gwas", "pvalue_eqtl", "effect_af_eqtl", "effect_af_gwas", "ZSCORE_gwas", "ZSCORE_eqtl"], index=False)

def launch_coloc(locus, window):
    # For quantitative traits there is no "s";
    # For case-control there may be no "beta" or "varbeta"
    if "s" not in locus.settings["gwas_experiments"][locus.gwas_file]:
        locus.settings["gwas_experiments"][locus.gwas_file]["s"] = 0

    if "N" not in locus.settings["eqtl_experiments"][locus.eqtl_file]:
        assert "eqtl_sample_sizes" in locus.settings
        eqtl_samples = pd.read_csv(locus.settings["eqtl_sample_sizes"], sep='\t', names=['eqtl_file','N'])
	eqtl_filename = locus.eqtl_file.split('/')[-1]
        N_eqtl = eqtl_samples.loc[eqtl_samples['eqtl_file']==eqtl_filename, "N"].iloc[0]
    else:
        N_eqtl = locus.settings["eqtl_experiments"][locus.eqtl_file]["N"]

    if "N" not in locus.settings["gwas_experiments"][locus.gwas_file]:
        assert "gwas_sample_sizes" in locus.settings
        gwas_samples = pd.read_csv(locus.settings["gwas_sample_sizes"], sep='\t', names=['gwas_file','N'])
        gwas_filename = locus.gwas_file.split('/')[-1]
        N_gwas = gwas_samples.loc[gwas_samples['gwas_file']==gwas_filename, "N"].iloc[0]
    else:
        N_gwas = locus.settings["gwas_experiments"][locus.gwas_file]["N"]


    # get "type"
    try:
        if "gwas_type" in locus.settings: # allow for global "type" value
            T = locus.settings["gwas_type"]
        else:
            T = locus.settings["gwas_experiments"][locus.gwas_file]["type"]
    except ValueError:
        print '"type" not specified for ' + locus.gwas_file
        raise

    # Launch an R script to run the coloc package.
    coloc_prob = subprocess.check_output("Rscript ./run_coloc.R {0}/coloc/{1}/{2}_{3}/{4}/{5}_level{6}.csv {7} {8} {9} {10}".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, \
        N_gwas, \
        locus.settings["gwas_experiments"][locus.gwas_file]["s"], \
        T, \
        N_eqtl), shell=True)
    # ^ Format of command line arguments: (N-gwas s-gwas type-gwas N-eqtl)

    coloc_prob_h0, coloc_prob_h1, coloc_prob_h2, coloc_prob_h3, coloc_prob_h4 = [float(x) for x in coloc_prob.strip().split()]

    num_sites = int(subprocess.check_output("wc -l {0}/coloc/{1}/{2}_{3}/{4}/{5}_level{6}.csv".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True).split()[0])-1

    # Write COLOC results to the desired file.
    with open("{0}/{1}_coloc_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, num_sites, coloc_prob_h0, coloc_prob_h1, coloc_prob_h2, coloc_prob_h3, coloc_prob_h4, locus.gwas_suffix))

    # Purge tmp files
    subprocess.call("rm -f {0}/coloc/{1}/{2}_{3}/{4}/{5}_level{6}.csv".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True)

    return coloc_prob_h4

# Get MAFs for data based on a reference VCF file
# Eliminate any data points for which MAFs can't be found
def get_mafs(locus, combined, window):

    ref = locus.settings["ref_genomes"][locus.settings["eqtl_experiments"][locus.eqtl_file]["ref"]]
    filename = ref["file"].format(locus.chrom)

    # determine presence of "chr" prefix in VCF
    if filename.endswith(".gz"):
        with gzip.open(filename, 'rb') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue
                if line.startswith("chr"):
                    chr_prefix = True
                    break
                else:
                    chr_prefix = False
                    break
    else:
        with open(filename, 'rb') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue
                if line.startswith("chr"):
                    chr_prefix = True
                    break
                else:
                    chr_prefix = False
                    break        

    # First, extract nearby variants using tabix
    #header = subprocess.check_output("zgrep -m1 -E '^#CHROM' {0} | cut -f1,2,3,4,5,8".format(filename), shell=True).strip().split()

    with gzip.open(filename, 'rb') as f:
        for line in f:
            if line.startswith("#CHROM"):
                header=[line.strip().split()[x] for x in [0,1,2,3,4,7]]
                break

    if chr_prefix:
        stream = StringIO(subprocess.check_output("tabix {8} {1}:{6}-{7} | cut -f1,2,3,4,5,8".format(locus.gwas_suffix, "chr" + str(locus.chrom), locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.pos-window, locus.pos+window, filename), shell=True))
    else:
        stream = StringIO(subprocess.check_output("tabix {8} {1}:{6}-{7} | cut -f1,2,3,4,5,8".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.pos-window, locus.pos+window, filename), shell=True))

    # For readability, load the header too
    # Load with pandas
    vcf = pd.read_csv(stream, sep="\t", names=header)

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

