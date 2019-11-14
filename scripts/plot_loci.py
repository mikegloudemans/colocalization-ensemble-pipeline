# Author: Mike Gloudemans
#
# Tools for plotting loci
#
# Produce a paired Manhattan plot and/or a colocalization
# scatterplot for the given locus. 

import subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab
import math
import numpy as np
import pandas as pd
import gzip
import sys
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO


def locus_compare(locus):

    data_copy = locus.data.copy()

    trait = locus.trait
    if trait == -1:
        trait = locus.gwas_suffix

    if "rsid_column" in locus.settings["gwas_experiments"][locus.gwas_file]:
        data_copy['rsid'] = data_copy[locus.settings["gwas_experiments"][locus.gwas_file]['rsid_column']]
    # If both GWAS and eQTL files have rsids, they'll have been renamed
    if "rsid_gwas" in list(data_copy.columns.values):
        data_copy['rsid'] = data_copy['rsid_gwas']

    if "rsid" not in list(data_copy.columns.values) or ("reobtain_rsids" in locus.settings and str(locus.settings["reobtain_rsids"]) == "True"):
        if "rsid_index_file" in locus.settings:

            # First, extract nearby variants using tabix
            stream = StringIO(subprocess.check_output("tabix {0} {1}:{2}-{3}".format(locus.settings["rsid_index_file"], str(data_copy['chr_gwas'][0]).replace("chr", ""), np.min(np.array(data_copy["snp_pos"])), np.max(np.array(data_copy["snp_pos"]))), shell=True))

            # For readability, load the header too
            # Load with pandas
            dbsnp = pd.read_csv(stream, sep="\t", comment="#", header=None).iloc[:,:3]
            dbsnp = dbsnp.rename({0: "chr_gwas", 1: "snp_pos", 2:"rsid"}, axis="columns")

            data_copy['chr_gwas'] = data_copy['chr_gwas'].astype(str).str.replace("chr", "", regex=False).astype(int)

            # May be unsafe, but I don't think it is because I don't think
            # we're using data_copy anywhere else
            data_copy = pd.merge(dbsnp, data_copy, left_on=["chr_gwas", "snp_pos"], right_on=["chr_gwas", "snp_pos"], suffixes = ["", "_original"])

            '''
            # Remove variants with position appearing multiple times
            dup_counts = {}
            for pos in vcf["POS"]:
                    dup_counts[pos] = dup_counts.get(pos, 0) + 1
            vcf["dup_counts"] = [dup_counts[pos] for pos in vcf['POS']]
            vcf = vcf[vcf["dup_counts"] == 1]
            '''

        else:
            return "No rsids found for plotting"

    # Throw away columns that don't have an rsid specified.
    plot_data = data_copy[~data_copy["rsid"].isna()]
 
    subprocess.call("mkdir -p {0}/plots/{4}/{5}/{1}_{2}/{3}".format(locus.basedir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix, trait), shell=True)
    gwas_out_file = "{0}/plots/{4}/{6}/{1}_{2}/{3}/{5}_gwas_locuscompare.png".format(locus.basedir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix, locus.gene, trait)
    eqtl_out_file = "{0}/plots/{4}/{6}/{1}_{2}/{3}/{5}_eqtl_locuscompare.png".format(locus.basedir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix, locus.gene, trait)

    # Assume that SNP rsid is already present, has been fetched earlier in program while
    # loading the GWAS
    subprocess.call(["mkdir", "-p", "{0}/locuscompare/{4}/{5}/{1}_{2}/{3}".format(locus.tmpdir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix, trait)])
    gwas_tmp = "{0}/locuscompare/{4}/{6}/{1}_{2}/{3}/{5}_gwas_lc_data.txt".format(locus.tmpdir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix, locus.gene, trait)
    eqtl_tmp = "{0}/locuscompare/{4}/{6}/{1}_{2}/{3}/{5}_eqtl_lc_data.txt".format(locus.tmpdir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix, locus.gene, trait)
 
    # Subset down to the region of interest, save this region
    vcf_file = "/mnt/lab_data/montgomery/shared/1KG/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz".format(locus.chrom)
    vcf_tmp = "{0}/locuscompare/{4}/{6}/{1}_{2}/{3}/{5}_vcf_tmp.vcf".format(locus.tmpdir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix, locus.gene, trait)
    with gzip.open(vcf_file) as f:
        line = f.readline()
        while line.startswith("##"):
            line = f.readline()
        with open(vcf_tmp, "w") as w:
            w.write(line)
    subprocess.check_call('tabix {3} {0}:{1}-{2} >> {4}'.format(locus.chrom, locus.pos - locus.settings["window"], locus.pos + locus.settings["window"], vcf_file, vcf_tmp), shell=True)
     
    gwas_data = data_copy.loc[:,["rsid", "pvalue_gwas"]]
    gwas_data.to_csv(gwas_tmp, header=["rsid", "pval"], index=False, sep="\t")
    eqtl_data = data_copy.loc[:,["rsid", "pvalue_eqtl"]]
    eqtl_data.to_csv(eqtl_tmp, header=["rsid", "pval"], index=False, sep="\t")

    gwas_data = gwas_data.reset_index()
    eqtl_data = eqtl_data.reset_index()

    gwas_lead = list(gwas_data["rsid"])[np.argmin(np.array(gwas_data["pvalue_gwas"]))]
    eqtl_lead = list(eqtl_data["rsid"])[np.argmin(np.array(eqtl_data["pvalue_eqtl"]))]

    # Call it once for top SNP in study 1, once for top SNP in study 2,
    # as reference SNP.
    # For now we'll just assume 1000 Genomes is the reference population
    subprocess.check_call(["Rscript", "./locuscompare.R", gwas_tmp, eqtl_tmp, gwas_out_file, locus.gwas_suffix, locus.eqtl_suffix, vcf_tmp, gwas_lead])
    subprocess.check_call(["Rscript", "./locuscompare.R", gwas_tmp, eqtl_tmp, eqtl_out_file, locus.gwas_suffix, locus.eqtl_suffix, vcf_tmp, eqtl_lead])

