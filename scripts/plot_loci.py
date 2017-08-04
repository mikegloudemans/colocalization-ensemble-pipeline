# Author: Mike Gloudemans
#
# Tools for plotting loci
#
# Produce a paired Manhattan plot and/or a colocalization
# scatterplot for the given locus. 

#import os
#import pandas as pd
#import sys
#import multiprocessing
#import gzip
#import collections

# Definitely need
import subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab
import math

def locus_zoom_plot(locus, clpp):
    gwas_chrom = locus.snp[0][3:]
    gwas_pos = locus.snp[1]
    locus.data = locus.data

    # For now
    current_level = 0

    subprocess.call("mkdir -p {0}/plots/{1}_{2}/{3}".format(locus.basedir, gwas_chrom, gwas_pos, locus.eqtl_suffix), shell=True)

    # Also create a LocusZoom-style plot showing the GWAS and eQTL signals next to one another.
    plt.figure(figsize=(20,10))
    plt.subplot(211)
    plt.scatter(locus.data['snp_pos'], [-1 * math.log10(p) for p in locus.data['pvalue_x']], c=locus.data['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)
    plt.plot((gwas_pos, gwas_pos), (-2, max([-1 * math.log10(p) for p in locus.data['pvalue_x']])), 'k--')
    plt.ylabel('GWAS -log p-value', fontsize=16)
    plt.title('{0} CLPP = {1}'.format(locus.gene, clpp), fontsize=24)
    plt.subplot(212)
    plt.scatter(locus.data['snp_pos'], [-1 * math.log10(p) for p in locus.data['pvalue_y']], c=locus.data['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)
    plt.plot((gwas_pos, gwas_pos), (-2, max([-1 * math.log10(p) for p in locus.data['pvalue_y']])), 'k--')
    plt.ylabel('eQTL -log p-value', fontsize=16)
    plt.xlabel('Position', fontsize=16)
    plt.savefig("{0}/plots/{1}_{2}/{3}/{4}_manhattan.png".format(locus.basedir, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene))
    plt.gcf().clear()
    plt.close()

def pvalue_plot(locus, clpp):
    gwas_chrom = locus.snp[0][3:]
    gwas_pos = locus.snp[1]

    # For now
    current_level = 0

    subprocess.call("mkdir -p {0}/plots/{1}_{2}/{3}".format(locus.basedir, gwas_chrom, gwas_pos, locus.eqtl_suffix), shell=True)

    plt.figure(figsize=(10,10))
    plt.scatter([-1 * math.log10(p) for p in locus.data['pvalue_x']], [-1 * math.log10(p) for p in locus.data['pvalue_y']], c=locus.data['snp_pos'], cmap=plt.cm.jet, edgecolor='', s=50)

    if max([-1 * math.log10(p) for p in locus.data['pvalue_x']] + [-1 * math.log10(p) for p in locus.data['pvalue_y']]) < 20:
            plt.axis([0, 20, 0, 20])
    else:
            plt.axis([0, max([20] + [-1 * math.log10(p) for p in locus.data['pvalue_x']] + [-1 * math.log10(p) for p in locus.data['pvalue_y']]), 0, max([20] + [-1 * math.log10(p) for p in locus.data['pvalue_x']] + [-1 * math.log10(p) for p in locus.data['pvalue_y']])])
    plt.xlabel('GWAS -log p-value', fontsize=16)
    plt.ylabel('eQTL -log p-value', fontsize=16)
    plt.title('{0} CLPP = {1}'.format(locus.gene, clpp), fontsize=24)
    plt.savefig("{0}/plots/{1}_{2}/{3}/{4}.png".format(locus.basedir, gwas_chrom, gwas_pos, locus.eqtl_suffix, locus.gene), shell=True)
    plt.gcf().clear()
    plt.close()

