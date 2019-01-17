# Author: Mike Gloudemans
#
# Run eCAVIAR algorithm for colocalization.
# Store full results in the output directory for this run.
# Return the CLPP score.
#

import subprocess
from scipy import stats
import math
import pandas as pd
import gzip
import os

def run_twas(locus, window=500000):

    subprocess.call("mkdir -p {4}/twas/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.tmpdir), shell=True)

    # Add rsids to the VCF files (may be unnecessary later if we change simulation format)
    # This is needed to match eQTL SNPs with GWAS SNPs
    vcf = pd.read_csv(locus.settings["eqtl_experiments"][locus.eqtl_file]["ref"], sep="\t", skiprows=1)
    gwas = pd.read_csv(locus.gwas_file, sep="\t")
    for i in range(vcf.shape[1]):
        assert gwas['variant_id'].iloc[i].startswith(vcf['ID'].iloc[i])
    vcf['ID'] = gwas['rsid']

    # We also need to clean out all the variants that don't appear in
    # our reference, so there won't be any weight assigned to them
    ref = pd.read_csv("/users/mgloud/software/TWAS-Fusion/fusion_twas-master/LDREF/nonbinary/1000G.EUR.{0}.map".format(locus.chrom), sep="\t", header=None)
    ref_list = set(ref.iloc[:,1])
    vcf = vcf[[v in ref_list for v in vcf['ID']]]

    vcf.to_csv('{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas.vcf'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), sep="\t", index=False)

    # Format genotypes data appropriately for TWAS to compute expression weights
    vcf_to_plink('{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas.vcf'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_plinked".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir))

    # Format phenotypes data appropriately for TWAS to compute expression weights
    add_phenos_to_fam("{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_plinked.fam".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_plinked.fam".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), locus.settings["eqtl_experiments"][locus.eqtl_file]["phenos"])
    
    # Compute expression weights
    # Several models are included...consider just deciding which is the best and using that. Or could
    # try and use all for evaluation?
    command = "Rscript /users/mgloud/software/TWAS-Fusion/fusion_twas-master/FUSION.compute_weights.R \
            --bfile {6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_plinked \
            --tmp {6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_tmp \
            --out {6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_weights \
            --models top1,lasso,enet \
            --hsq_p 1 \
            --verbose 0".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)
    subprocess.check_call(command, shell=True)
  
    # Some genes don't make it through the above part, presumably because there's no
    # strong eQTL. If we can't build a model for the gene, then I guess we just have to skip it for now.
    # Later we might set the --hsq_p flag above to just let everything pass, and see if it does worse that way
    if not os.path.isfile("{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_weights.wgt.RDat".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)):
        # If the previous step failed to produce weights, there's not much we can do after that
        return 0

    # Format GWAS correctly for TWAS to run
    gwas = pd.read_csv(locus.gwas_file, sep="\t")
    gwas = gwas[['rsid', 'alt', 'ref', 'zscore']]
    gwas.columns = ['SNP', 'A1', 'A2', 'Z']
    gwas.to_csv("{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), sep="\t", index=False)
    
    # Make POS file
    with open("{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_weights.pos".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "w") as w:
        w.write("WGT\tID\tCHR\tP0\tP1\n")
        w.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format("{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_weights.wgt.RDat".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "nameless_gene", locus.chrom, locus.pos, locus.pos))

    # Run TWAS with permutations
    command = "Rscript /users/mgloud/software/TWAS-Fusion/fusion_twas-master/FUSION.assoc_test.R \
            --sumstats {6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas \
            --weights {6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_weights.pos \
            --weights_dir / \
            --ref_ld_chr /users/mgloud/software/TWAS-Fusion/fusion_twas-master/LDREF/1000G.EUR. \
            --chr {1} \
            --out {6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_results.dat \
            --perm 1000".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)
    subprocess.check_call(command, shell=True)

    # Parse results
    print "{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_results.dat".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)
    with open("{6}/twas/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_twas_results.dat".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)) as f:
        f.readline()
        data = f.readline().strip().split()
        twas_p = float(data[18])
        #twas_p_perm = float(data[21]) # Figure out what the deal with this is before trying to report it and crashing the script
        twas_p_perm = "NA"
        snps_tested = int(data[12])

    print twas_p, twas_p_perm, snps_tested

    # Add results to the desired file
    with open("{0}/{1}_twas_clpp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, snps_tested, -1*math.log10(twas_p), twas_p_perm))
    
    return twas_p


def vcf_to_plink(vcf_file, plink_file):
    # Use PLINK to generate bim bam fam files
    command = '''/srv/persistent/bliu2/tools/plink_1.90_beta3_linux_x86_64/plink --vcf {0} --keep-allele-order --make-bed --double-id --out {1} > /dev/null'''.format(vcf_file, plink_file)
    subprocess.check_call(command, shell=True)

def add_phenos_to_fam(old_fam, new_fam, pheno_file):

    with gzip.open(pheno_file) as f:
        phenos = pd.read_csv(f, sep="\t")
        phenos = list(phenos.iloc[0, 6:])

    table = pd.read_csv(old_fam, sep = " ", header=None)
    table.iloc[:, 5] = phenos
    table.to_csv(new_fam, sep=" ", header=None, index=False)
