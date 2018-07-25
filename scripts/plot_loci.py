# Author: Mike Gloudemans
#
# Tools for plotting loci
#
# Produce a paired Manhattan plot and/or a colocalization
# scatterplot for the given locus. 

# Definitely need
import subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab
import math
import numpy as np
import pandas as pd

def locus_zoom_plot(locus, clpp):

    subprocess.call("mkdir -p {0}/plots/{4}/{1}_{2}/{3}".format(locus.basedir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix), shell=True)

    # Also create a LocusZoom-style plot showing the GWAS and eQTL signals next to one another.
    plt.figure(figsize=(20,10))
    plt.subplot(211)

    plt.scatter(locus.data['snp_pos'], [-1 * math.log10(p) for p in locus.data['pvalue_gwas']], c=locus.data['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)
    plt.plot((locus.pos, locus.pos), (-0.1*max([-1 * math.log10(p) for p in locus.data['pvalue_gwas']]), max([-1 * math.log10(p) for p in locus.data['pvalue_gwas']])), 'k--')
    plt.ylabel('GWAS -log p-value', fontsize=16)
    plt.title('{0} CLPP = {1}'.format(locus.gene, clpp), fontsize=24)
    plt.subplot(212)
    plt.scatter(locus.data['snp_pos'], [-1 * math.log10(p) for p in locus.data['pvalue_eqtl']], c=locus.data['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)
    plt.plot((locus.pos, locus.pos), (-0.1*max([-1 * math.log10(p) for p in locus.data['pvalue_eqtl']]), max([-1 * math.log10(p) for p in locus.data['pvalue_eqtl']])), 'k--')
    plt.ylabel('eQTL -log p-value', fontsize=16)
    plt.xlabel('Position', fontsize=16)
    plt.savefig("{0}/plots/{5}/{1}_{2}/{3}/{4}_manhattan.png".format(locus.basedir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.gwas_suffix))
    plt.gcf().clear()
    plt.close()

def pvalue_plot(locus, clpp):

    subprocess.call("mkdir -p {0}/plots/{4}/{1}_{2}/{3}".format(locus.basedir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix), shell=True)

    plt.figure(figsize=(10,10))
    plt.scatter([-1 * math.log10(p) for p in locus.data['pvalue_gwas']], [-1 * math.log10(p) for p in locus.data['pvalue_eqtl']], c=locus.data['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)

    if max([-1 * math.log10(p) for p in locus.data['pvalue_gwas']] + [-1 * math.log10(p) for p in locus.data['pvalue_eqtl']]) < 20:
            plt.axis([0, 20, 0, 20])
    else:
            plt.axis([0, max([20] + [-1 * math.log10(p) for p in locus.data['pvalue_gwas']] + [-1 * math.log10(p) for p in locus.data['pvalue_eqtl']]), 0, max([20] + [-1 * math.log10(p) for p in locus.data['pvalue_gwas']] + [-1 * math.log10(p) for p in locus.data['pvalue_eqtl']])])
    plt.xlabel('GWAS -log p-value', fontsize=16)
    plt.ylabel('eQTL -log p-value', fontsize=16)
    plt.title('{0} CLPP = {1}'.format(locus.gene, clpp), fontsize=24)
    plt.savefig("{0}/plots/{5}/{1}_{2}/{3}/{4}.png".format(locus.basedir, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.gwas_suffix), shell=True)
    plt.gcf().clear()
    plt.close()

def locus_compare(locus):

    trait = locus.trait
    if trait == -1:
        trait = locus.gwas_suffix

    if "rsid_column" in locus.settings["gwas_experiments"][locus.gwas_file]:
        locus.data['rsid'] = locus.data[locus.settings["gwas_experiments"][locus.gwas_file]['rsid_column']]
    # If both GWAS and eQTL files have rsids, they'll have been renamed
    if "rsid_gwas" in list(locus.data.columns.values):
        locus.data['rsid'] = locus.data['rsid_gwas']

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
    subprocess.check_call("zcat {0} | tail -n +253 | head -n 1 > {1}".format(vcf_file, vcf_tmp), shell=True)
    subprocess.check_call('tabix {3} {0}:{1}-{2} >> {4}'.format(locus.chrom, locus.pos - locus.settings["window"], locus.pos + locus.settings["window"], vcf_file, vcf_tmp), shell=True)
     
    gwas_data = locus.data.loc[:,["rsid", "pvalue_gwas"]]
    gwas_data.to_csv(gwas_tmp, header=["rsid", "pval"], index=False, sep="\t")
    eqtl_data = locus.data.loc[:,["rsid", "pvalue_eqtl"]]
    eqtl_data.to_csv(eqtl_tmp, header=["rsid", "pval"], index=False, sep="\t")

    gwas_data = gwas_data.reset_index()
    eqtl_data = eqtl_data.reset_index()

    gwas_lead = list(gwas_data["rsid"])[np.argmin(gwas_data["pvalue_gwas"])]
    eqtl_lead = list(eqtl_data["rsid"])[np.argmin(eqtl_data["pvalue_eqtl"])]

    # Call it once for top SNP in study 1, once for top SNP in study 2,
    # as reference SNP.
    # For now we'll just assume 1000 Genomes is the reference population
    print "Rscript", "/users/mgloud/projects/brain_gwas/scripts/locuscompare.R", gwas_tmp, eqtl_tmp, gwas_out_file, locus.gwas_suffix, locus.eqtl_suffix, vcf_tmp, gwas_lead
    subprocess.check_call(["Rscript", "/users/mgloud/projects/brain_gwas/scripts/locuscompare.R", gwas_tmp, eqtl_tmp, gwas_out_file, locus.gwas_suffix, locus.eqtl_suffix, vcf_tmp, gwas_lead])
    subprocess.check_call(["Rscript", "/users/mgloud/projects/brain_gwas/scripts/locuscompare.R", gwas_tmp, eqtl_tmp, eqtl_out_file, locus.gwas_suffix, locus.eqtl_suffix, vcf_tmp, eqtl_lead])

