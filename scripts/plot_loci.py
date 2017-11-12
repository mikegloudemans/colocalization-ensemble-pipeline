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

def locus_zoom_plot(locus, clpp):

    # For now
    current_level = 0

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
    # For now
    current_level = 0

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

