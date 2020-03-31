# Author: Mike Gloudemans
#
# Run FINEMAP algorithm for colocalization.
# Store full results in the output directory for this run.
# Return the CLPP score.
#

import sys
import subprocess
from scipy import stats
from shutil import copyfile
if sys.version_info[0] < 3: 
   from StringIO import StringIO
else:
   from io import StringIO
import pandas as pd
import numpy as np
import math
from operator import mul
import gzip

# TODO: Make it so that different traits are written in different temporary files
# Otherwise there may be concurrency bugs if running in parallel across loci

def run_finemap(locus, window=500000):

    pf = prep_finemap(locus, window)

    # TODO: Make it write to an error or a "skipped" file instead of failing silently.
    if isinstance(pf, basestring):
        print pf
        return pf
    return launch_finemap(locus, window, pf)


# Separates the preparation of data for finemap from the actual 
# launching of finemap. This is done to avoid duplicating code,
# since eCAVIAR and FINEMAP use very similar setups.
def prep_finemap(locus, window):

    locus.conditional_level = 0   # Currently not used at all, but we may re-add this functionality later.
    combined = locus.data.copy()

    # Make required directories for FINEMAP eCAVIAR analysis
    subprocess.call("mkdir -p {0}/vcftools/{1}/{2}_{3}/{4}".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
    subprocess.call("mkdir -p {0}/plink/{1}/{2}_{3}/{4}".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
    subprocess.call("mkdir -p {0}/ecaviar/{1}/{2}_{3}/{4}".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
    subprocess.call("mkdir -p {0}/finemap".format(locus.basedir), shell=True)
    
    # Get VCF file paths
    eqtl_ref = locus.settings["ref_genomes"][locus.settings["eqtl_experiments"][locus.eqtl_file]["ref"]]
    gwas_ref = locus.settings["ref_genomes"][locus.settings["gwas_experiments"][locus.gwas_file]["ref"]]
    eqtl_vcf = eqtl_ref["file"].format(locus.chrom)
    gwas_vcf = gwas_ref["file"].format(locus.chrom)

    # Two different cases depending on whether GWAS and eQTL
    # are using same reference genome.
    if eqtl_vcf == gwas_vcf:
        # Get and filter the single VCF.

        vcf, combined = load_and_filter_variants(eqtl_vcf, locus, combined, eqtl_ref, window, ["eqtl", "gwas"])
        assert vcf.shape[0] == combined.shape[0]

        # Run PLINK on just one VCF.
        removal_list = compute_ld(vcf, locus, "eqtl")
        if isinstance(removal_list, basestring):
            return removal_list
        subprocess.check_call("cp {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.fixed.ld {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.fixed.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), shell=True)

        # Remove indices that produced NaNs in the LD computations
        removal_list = list(set(removal_list))
        combined = combined.drop(combined.index[removal_list])
    else:
        # Get and filter both VCFs.
        evcf, combined = load_and_filter_variants(eqtl_vcf, locus, combined, eqtl_ref, window, ["eqtl"])
        gvcf, combined = load_and_filter_variants(gwas_vcf, locus, combined, gwas_ref, window, ["gwas"])

        # Subset to overlapping SNPs
        evcf, gvcf, combined = intersect_reference_vcfs(evcf, gvcf, combined)
        assert evcf.shape[0] == combined.shape[0]
        assert gvcf.shape[0] == combined.shape[0]

        # Run PLINK on both VCFs.
        while True:
            removal_list = compute_ld(evcf, locus, "eqtl")
            if isinstance(removal_list, basestring):
                return removal_list
            
            extension_list = compute_ld(gvcf, locus, "gwas")
            if isinstance(extension_list, basestring):
                return extension_list
 
            removal_list.extend(extension_list)

            # Continue until no more NaNs
            if len(removal_list) == 0:
                break

            removal_list = list(set(removal_list))

            combined = combined.drop(combined.index[removal_list])
            evcf = evcf.drop(evcf.index[removal_list])
            gvcf = gvcf.drop(gvcf.index[removal_list])


    # Check to see whether we still even have a signficant GWAS variant and a significant eQTL variant.
    if "screening_thresholds" in locus.settings and "gwas" in locus.settings["screening_thresholds"]:
        if min(combined["pvalue_gwas"]) > locus.settings["screening_thresholds"]["gwas"]:
            return "Fail: No GWAS SNPs pass thresholds after filtering."
    if "screening_thresholds" in locus.settings and "eqtl" in locus.settings["screening_thresholds"]:
        if min(combined["pvalue_eqtl"]) > locus.settings["screening_thresholds"]["eqtl"]:
            return "Fail: No eQTL SNPs pass thresholds after filtering."
    
    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "w") as w:
        snps = combined[['snp_pos', 'ZSCORE_eqtl']]
        snps.to_csv(w, index=False, header=False, sep=" ")

    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "w") as w:
        snps = combined[['snp_pos', 'ZSCORE_gwas']]
        snps.to_csv(w, index=False, header=False, sep=" ")
    return (min(combined["pvalue_gwas"]), min(combined["pvalue_eqtl"]))

# This function contains the code that's specific to FINEMAP,
# not shared with eCAVIAR.
def launch_finemap(locus, window, top_hits):

    # Load sample sizes
    eqtl_n = locus.settings["ref_genomes"][locus.settings["eqtl_experiments"][locus.eqtl_file]["ref"]]["N"]
    gwas_n = locus.settings["ref_genomes"][locus.settings["gwas_experiments"][locus.gwas_file]["ref"]]["N"]

    # Write config file for finemap
    subprocess.check_call('echo "z;ld;snp;config;n-ind" > {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), shell=True)
    subprocess.check_call('echo "{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.fixed.ld;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.config;{6}" >> {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, gwas_n, locus.tmpdir), shell=True)
    subprocess.check_call('echo "{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.fixed.ld;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.config;{6}" >> {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, eqtl_n, locus.tmpdir), shell=True)
   
    max_causal = 1
    if "max_causal" in locus.settings["methods"]["finemap"]:
        max_causal = int(locus.settings["methods"]["finemap"]["max_causal"])

    # Run FINEMAP
    if locus.settings["debug"] == True:
        subprocess.check_call('finemap --sss --in-files {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in --n-causal-snps {7} --n-iterations 1000000 --n-conv-sss 1000'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir, max_causal), shell=True)
    else:
        subprocess.check_call('finemap --sss --in-files {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in --n-causal-snps {7} --n-iterations 1000000 --n-conv-sss 1000 > /dev/null'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir, max_causal), shell=True)
    
    # Parse FINEMAP results to compute CLPP score
    gwas_probs = []
    eqtl_probs = []
    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)) as f:
            f.readline()
            for line in f:
                    data = line.strip().split()
                    gwas_probs.append((int(data[1]), float(data[2])))
    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)) as f:
            f.readline()
            for line in f:
                    data = line.strip().split()
                    eqtl_probs.append((int(data[1]), float(data[2])))
    gwas_probs = sorted(gwas_probs)
    eqtl_probs = sorted(eqtl_probs)

    assert len(gwas_probs) == len(eqtl_probs)
    for i in range(len(gwas_probs)):
            assert gwas_probs[i][0] == eqtl_probs[i][0]

    #finemap_clpp = sum([gwas_probs[i][1] * eqtl_probs[i][1] for i in range(len(gwas_probs))])
    
    finemap_clpp = 1 - reduce(mul, [1-(gwas_probs[i][1]*eqtl_probs[i][1]) for i in range(len(gwas_probs))])

    ld_file = "{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.fixed.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)
    finemap_clpp_mod = get_clpp_mod(gwas_probs, eqtl_probs, ld_file)

    # Write header of output file for FINEMAP
    trait_suffix = locus.trait.split("/")[-1].replace(".", "_")

    if finemap_clpp_mod > 0.3 or ("save_all_finemap" in locus.settings["methods"]["finemap"] and locus.settings["methods"]["finemap"]["save_all_finemap"] == "True"):
    #if finemap_clpp_mod > 0:
        copyfile("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "{0}/finemap/{1}_{2}_{3}_{4}_{5}_{6}_{7}_finemap_gwas.snp".format(locus.basedir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, trait_suffix))
        copyfile("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "{0}/finemap/{1}_{2}_{3}_{4}_{5}_{6}_{7}_finemap_eqtl.snp".format(locus.basedir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, trait_suffix))

    # Write FINEMAP results to the desired file
    # Note: appending will always work, since results always go to a different directory for each run.
    # TODO: Write headers for this file
    with open("{0}/{1}_finemap_clpp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
                a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, len(gwas_probs), finemap_clpp, -1*math.log10(top_hits[0]), -1*math.log10(top_hits[1]), locus.gwas_suffix, finemap_clpp_mod))

    return finemap_clpp, finemap_clpp_mod

# Function purge_tmp_files
# Remove temporary files created during this run of eCAVIAR,
# to free up space.
def purge_tmp_files(locus):
    
    # Why not just purge everything from the entire directory of the GWAS that we're working on?
    # Answer: it's possible that other jobs are still using those files.
    # Only purge the specific files created in this run.

    # (To avoid conflicts here, just make sure no two processes are examining exact same SNP/gene
    # pair at the same time.)

    subprocess.call("rm -rf {4}/vcftools/{0}/{1}_{2}/{3}/{5}*".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.tmpdir, locus.gene), shell=True)
    subprocess.call("rm -rf {4}/ecaviar/{0}/{1}_{2}/{3}/{5}*".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.tmpdir, locus.gene), shell=True)
    subprocess.call("rm -rf {4}/plink/{0}/{1}_{2}/{3}/{5}*".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.tmpdir, locus.gene), shell=True)

def load_and_filter_variants(filename, locus, combined, ref, window, ref_types):
    
    # First, extract nearby variants using tabix
    with gzip.open(filename) as f:
        line = f.readline()
        while line.startswith("##"):
            line = f.readline()
        header = line.strip().split()
      
    # get the proper chromosome prefix
    if "chr_prefix" in ref:
        chrom_prefix = ref["chr_prefix"]
    else:
        # pull it from the first non-header line of the VCF 
        with gzip.open(filename) as f:
            while line.startswith("#"):
                line = f.readline()
                  
        l = line.strip().split()
        chrom = l[0]
        if chrom.startswith('chr'):
            chrom_prefix = 'chr'
        elif chrom[0].isdigit():
            chrom_prefix = '' 
   
    stream = StringIO(subprocess.check_output("tabix {8} {1}:{6}-{7}".format(locus.gwas_suffix, chrom_prefix + str(locus.chrom), locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.pos-window, locus.pos+window, filename), shell=True))
   
    #if "chr_prefix" in ref and ref["chr_prefix"] == "chr":
    #    stream = StringIO(subprocess.check_output("tabix {8} {1}:{6}-{7}".format(locus.gwas_suffix, "chr" + str(locus.chrom), locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.pos-window, locus.pos+window, filename), shell=True))
    #else:
    #    stream = StringIO(subprocess.check_output("tabix {8} {1}:{6}-{7}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.pos-window, locus.pos+window, filename), shell=True))

    # For readability, load the header too
    # Load with pandas
    vcf = pd.read_csv(stream, sep="\t", names=header)
    print 'Number of SNPs in VCF: {}'.format(vcf.shape)
    # Remove variants not in the GWAS table
    vcf["POS"] = (vcf["POS"]).astype(int)
    vcf = vcf[vcf["POS"].isin(list(combined["snp_pos"]))]
    print 'Number of SNPs in VCF after filtering on combined["snp_pos"]: {}'.format(vcf.shape)

    # Remove variants with position appearing multiple times
    dup_counts = {}
    for pos in vcf["POS"]:
            dup_counts[pos] = dup_counts.get(pos, 0) + 1
    vcf["dup_counts"] = [dup_counts[pos] for pos in vcf['POS']]
    vcf = vcf[vcf["dup_counts"] == 1]
    print 'Number of SNPs in VCF after removing dup_counts: {}'.format(vcf.shape)

    # Remove multiallelic variants with only one entry in VCF
    l = lambda x: "," not in x
    vcf = vcf[vcf["REF"].apply(l) & vcf["ALT"].apply(l)]
    print 'Number of SNPs in VCF after removing multiallelic variants: {}'.format(vcf.shape)

    # Remove monoallelic variants.
    # Allele frequency might be input as counts or as percentages,
    # so handle this.
    if "filtered_by_af" not in ref or ref["filtered_by_af"] != "True":
        if "af_attribute" in ref:
            af_id = ref["af_attribute"]
            def fn(x):
                info = [s for s in x.split(";") if s.startswith(af_id + "=")][0]
                af = float(info.split("=")[1])
                return af > 0.01 and 1-af > 0.01 
            vcf = vcf[vcf["INFO"].apply(fn)]
            print 'Number of SNPs in VCF after filtering on AF: {}'.format(vcf.shape)
        else:
            ac_id = ref["ac_attribute"]
            an = 2*ref["N"] # Assume 2 alleles per person
            def fn(x):
                info = [s for s in x.split(";") if s.startswith(ac_id + "=")][0]
                ac = float(info.split("=")[1])
                af = ac*1.0/an
                return af > 0.01 and 1-af > 0.01 
            vcf = vcf[vcf["INFO"].apply(fn)]

    if "POS" in list(combined.columns.values):
        combined = combined.drop(columns=['POS'])

    # Remove variants where alt/ref don't match between GWAS/eQTL and VCF
    # Flipped is okay. A/C and C/A are fine, A/C and A/G not fine.
    # TODO: Verify on an example case that this filtering is working correctly
    print 'Number of SNPs in VCF before merging: {}'.format(vcf.shape)
    print 'Number of SNPs in "combined" before merging: {}'.format(combined.shape)
    merged = pd.merge(combined, vcf, left_on="snp_pos", right_on="POS")

    # TODO: Enforce new standard: effect measurements are always with respect to ALT status.
    keep_indices = \
            (((merged['ref_gwas'] == merged['REF']) & (merged['alt_gwas'] == merged['ALT'])) | \
            ((merged['alt_gwas'] == merged['REF']) & (merged['ref_gwas'] == merged['ALT']))) & \
            (((merged['ref_eqtl'] == merged['REF']) & (merged['alt_eqtl'] == merged['ALT'])) | \
            ((merged['alt_eqtl'] == merged['REF']) & (merged['ref_eqtl'] == merged['ALT'])))

    # Now, reverse the z-score for any SNPs whose ref/alt direction doesn't match the direction
    # in the reference genome.
    assert "gwas" in ref_types or "eqtl" in ref_types
    if "gwas" in ref_types:
        merged['reverse_gwas'] = (merged['alt_gwas'] == merged['REF']) & (merged['ref_gwas'] == merged['ALT'])
        merged['ZSCORE_gwas'] = merged['ZSCORE_gwas'] * (1 - (merged['reverse_gwas'] * 2))
    if "eqtl" in ref_types:
        merged['reverse_eqtl'] = (merged['alt_eqtl'] == merged['REF']) & (merged['ref_eqtl'] == merged['ALT'])
        merged['ZSCORE_eqtl'] = merged['ZSCORE_eqtl'] * (1 - (merged['reverse_eqtl'] * 2))

    keep = merged['POS'][keep_indices]
    vcf = vcf[vcf['POS'].isin(list(keep))]
    
    # Subset SNPs down to SNPs present in the reference VCF.
    combined = combined[combined['snp_pos'].isin(list(vcf["POS"]))]

    # Return list as DataFrame.
    return (vcf, combined)

def intersect_reference_vcfs(ref1, ref2, combined):

    # Subset each reference VCF down to only variants present in the other VCF.
    ref1 = ref1[ref1["POS"].isin(list(ref2["POS"]))]
    ref2 = ref2[ref2["POS"].isin(list(ref1["POS"]))]

    # Subset SNP list to SNPs that still remain in both VCFs.
    combined = combined[combined['snp_pos'].isin(list(ref2["POS"]))]

    return (ref1, ref2, combined)

# Run PLINK on the locus of interest
def compute_ld(input_vcf, locus, data_type):

    # We don't want to modify the input VCF within this function
    vcf = input_vcf.copy()

    # Repeatedly compute LD until we've eliminated all NaNs.
    removal_list = []
    while True:
        if vcf.shape[0] == 0:
            return "Fail: All SNPs have been eliminated through VCF filtering."

        # Write VCF to tmp file
        vcf.to_csv('{7}/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.{6}.vcf'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir), sep="\t", index=False, header=True)

        # Use PLINK to generate bim bam fam files
        command = '''plink --vcf {7}/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.{6}.vcf --keep-allele-order --make-bed --double-id --out {7}/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}_plinked > /dev/null'''.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir)
        subprocess.check_call(command, shell=True)

        # Use PLINK to generate LD score
        command = '''plink -bfile {7}/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}_plinked --r square --out {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6} > /dev/null'''.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir)
        subprocess.check_call(command, shell=True)

        # See if nans remain. If so, remove the offending lines.
        try:
            lines = [int(n.split(":")[0])-1 for n in subprocess.check_output("grep -n nan {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir), shell=True).strip().split("\n")]
        except:
            break


        # Save IDs of items being removed
        removal_list.extend(list(vcf.iloc[lines]['ID']))
        
        # Remove desired rows (variants)
        vcf = vcf.drop(vcf.index[lines])

    # Get indices of items to remove in original list
    removal_list = [list(input_vcf['ID']).index(id) for id in removal_list]

    # Replace tabs with spaces because FINEMAP requires this.
    subprocess.check_call("cat {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.ld | sed s/\\\\t/\\ /g > {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.fixed.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir), shell=True)

    return removal_list

# In this function, gwas and eqtl probs should be pre-sorted
def get_clpp_mod(gwas_probs, eqtl_probs, ld_file):
    
    ld = []
    with open(ld_file) as f:
        for line in f:
            ld.append([float(f) for f in line.strip().split()]) 


    for i in range(len(gwas_probs)):
        assert gwas_probs[i][0] == eqtl_probs[i][0]

    # Get modifeid CLPP score
    clpp_mod = 0
    for i in range(len(gwas_probs)):
        for j in range(len(eqtl_probs)):
            snp_ld = ld[i][j]**2
            clpp_mod += gwas_probs[i][1] * eqtl_probs[j][1] * snp_ld
    
    return clpp_mod
