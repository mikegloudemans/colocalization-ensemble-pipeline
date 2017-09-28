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

def run_finemap(locus, window=500000):

    prep_finemap(locus, window)
    return launch_finemap(locus, window)


# Separates the preparation of data for finemap from the actual 
# launching of finemap. This is done to avoid duplicating code,
# since eCAVIAR and FINEMAP use very similar setups.
def prep_finemap(locus, window):

    locus.conditional_level = 0   # Currently not used at all, but we may re-add this functionality later.
    combined = locus.data.copy()

    # Make required directories for FINEMAP eCAVIAR analysis
    subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
    subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
    subprocess.call("mkdir -p /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
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
        vcf, combined = load_and_filter_variants(eqtl_vcf, locus, combined, eqtl_ref, window)
        assert vcf.shape[0] == combined.shape[0]

        # Run PLINK on just one VCF.
        compute_ld(vcf, locus, "eqtl")
        subprocess.check_call("cp /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.fixed.ld /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.fixed.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True)


    else:
        # Get and filter both VCFs.
        evcf, combined = load_and_filter_variants(eqtl_vcf, locus, combined, eqtl_ref, window)
        gvcf, combined = load_and_filter_variants(gwas_vcf, locus, combined, gwas_ref, window)

        # Subset to overlapping SNPs
        evcf, gvcf, combined = intersect_reference_vcfs(evcf, gvcf, combined)
        assert evcf.shape[0] == combined.shape[0]
        assert gvcf.shape[0] == combined.shape[0]

        # Run PLINK on both VCFs.
        compute_ld(evcf, locus, "eqtl")
        compute_ld(gvcf, locus, "gwas")
    
    with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), "w") as w:
        snps = combined[['snp_pos', 'ZSCORE_eqtl']]
        snps.to_csv(w, index=False, header=False, sep=" ")

    with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), "w") as w:
        snps = combined[['snp_pos', 'ZSCORE_gwas']]	
        snps.to_csv(w, index=False, header=False, sep=" ")
    


# This function contains the code that's specific to FINEMAP,
# not shared with eCAVIAR.
def launch_finemap(locus, window):

    # Load sample sizes
    eqtl_n = locus.settings["ref_genomes"][locus.settings["eqtl_experiments"][locus.eqtl_file]["ref"]]["N"]
    gwas_n = locus.settings["ref_genomes"][locus.settings["gwas_experiments"][locus.gwas_file]["ref"]]["N"]

    # Write config file for finemap
    subprocess.check_call('echo "z;ld;snp;config;n-ind" > /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True)
    subprocess.check_call('echo "/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.fixed.ld;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.config;{6}" >> /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, gwas_n), shell=True)
    subprocess.check_call('echo "/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.fixed.ld;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp;/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.config;{6}" >> /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, eqtl_n), shell=True)
    
    # Run FINEMAP
    subprocess.check_call('finemap --sss --in-files /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in --n-causal-max 1 --n-iterations 1000000 --n-convergence 50000 2> /dev/null'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), shell=True)

    # Parse FINEMAP results to compute CLPP score
    gwas_probs = []
    eqtl_probs = []
    with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level)) as f:
            f.readline()
            for line in f:
                    data = line.strip().split()
                    gwas_probs.append((int(data[1]), float(data[2])))
    with open("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level)) as f:
            f.readline()
            for line in f:
                    data = line.strip().split()
                    eqtl_probs.append((int(data[1]), float(data[2])))
    gwas_probs = sorted(gwas_probs)
    eqtl_probs = sorted(eqtl_probs)



    assert len(gwas_probs) == len(eqtl_probs)
    for i in range(len(gwas_probs)):
            assert gwas_probs[i][0] == eqtl_probs[i][0]

    finemap_clpp = sum([gwas_probs[i][1] * eqtl_probs[i][1] for i in range(len(gwas_probs))])

    if finemap_clpp > 0.001:
        copyfile("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), "{0}/finemap/{1}_{2}_{3}_{4}_{5}_{6}_finemap_gwas.snp".format(locus.basedir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level))
        copyfile("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level), "{0}/finemap/{1}_{2}_{3}_{4}_{5}_{6}_finemap_eqtl.snp".format(locus.basedir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level))

    # Write FINEMAP results to the desired file
    # Note: appending will always work, since results always go to a different directory for each run.
    with open("{0}/{1}_finemap_clpp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
            a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, len(gwas_probs), finemap_clpp))

    return finemap_clpp

# Function purge_tmp_files
# Remove temporary files created during this run of eCAVIAR,
# to free up space.
def purge_tmp_files(locus):
   
    # Why not just purge everything from the entire directory of the GWAS that we're working on?
    # Answer: it's possible that other jobs are still using those files.
    # Only purge the specific files created in this run.

    # (To avoid conflicts here, make sure each gwas/eqtl/snp triplet is run in its own process;
    # doubling up would cause problems).

    subprocess.call("rm -r /users/mgloud/projects/brain_gwas/tmp/vcftools/{0}/{1}_{2}/{3} 2> /dev/null".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene), shell=True)
    subprocess.call("rm -r /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3} 2> /dev/null".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene), shell=True)
    subprocess.call("rm -r /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3} 2> /dev/null".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene), shell=True)

def load_and_filter_variants(filename, locus, combined, ref, window):
    
    # TODO: Spot check all of these tests to ensure they're working as desired.

    # First, extract nearby variants using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 500 | grep \\#CHROM".format(filename), shell=True).strip().split()
    stream = StringIO(subprocess.check_output("tabix {8} {1}:{6}-{7}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.pos-window, locus.pos+window, filename), shell=True))

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
            return af > 0.01 and 1-af > 0.01 
        vcf = vcf[vcf["INFO"].apply(fn)]
    else:
        ac_id = ref["ac_attribute"]
        an = 2*ref["N"] # Assume 2 alleles per person
        def fn(x):
            info = [s for s in x.split(";") if s.startswith(ac_id + "=")][0]
            ac = float(info.split("=")[1])
            af = ac*1.0/an
            return af > 0.01 and 1-af > 0.01 
        vcf = vcf[vcf["INFO"].apply(fn)]

    # Remove variants where alt/ref don't match between GWAS and VCF
    # Flipped is okay. A/C and C/A are fine, A/C and A/G not fine.
    # TODO: Verify on an example case that this filtering is working correctly.
    merged = pd.merge(combined, vcf, left_on="snp_pos", right_on="POS")
    keep_indices = ((merged['a1'] == merged['REF']) & (merged['a2'] == merged['ALT'])) | \
            ((merged['a2'] == merged['REF']) & (merged['a1'] == merged['ALT']))
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
def compute_ld(vcf, locus, data_type):


    # Write VCF to tmp file
    vcf.to_csv('/users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.{6}.vcf'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type), sep="\t", index=False, header=True)

    # Use PLINK to generate bim bam fam files
    command = '''/srv/persistent/bliu2/tools/plink_1.90_beta3_linux_x86_64/plink --vcf /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.{6}.vcf --keep-allele-order --make-bed --double-id --out /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}_plinked 2> /dev/null'''.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type)
    subprocess.check_call(command, shell=True)

    # Use PLINK to generate LD score
    command = '''/srv/persistent/bliu2/tools/plink_1.90_beta3_linux_x86_64/plink -bfile /users/mgloud/projects/brain_gwas/tmp/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}_plinked --r square --out /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6} 2> /dev/null'''.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type)
    subprocess.check_call(command, shell=True)

    # Monitor nan status. Goal is eventually to just get rid of this, either
    # because I have no nans or because I completely understand them
    try:
        subprocess.check_call("grep nan /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type), shell=True)
        # Log scenarios with nans to a file.
        with open("{0}/nans.txt".format(locus.basedir), "a") as w:
            w.write("/users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.ld\n".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level))
    except:
        pass

    # Replace tabs with spaces because FINEMAP requires this.
    # Fix LD-score by replacing nan values with 0.
    # TODO: If nans are inevitable, maybe replace them with imputed values instead of 0s.
    # As it is currently, it could affect results in certain rare cases.
    # Figure out why these are nanning in the first place
    subprocess.check_call("sed s/nan/0/g /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.ld | sed s/\\\\t/\\ /g > /users/mgloud/projects/brain_gwas/tmp/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.fixed.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type), shell=True)
