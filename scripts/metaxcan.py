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

def run_metaxcan(locus, window=500000):

    ###########################################################
    # Using eQTLs and phenotypes, generate the weights database
    ###########################################################

    # Transform gene expression data into the format required for PredictDBPipeline
    # This is easy; it's just a matrix

    # Make a basic GENCODE file with just the one annotation we're interested in
    with open("{0}/gencode.txt".format(tmpdir), "w") as w:
        w.write('''{0}\tHAVANA\tgene\t{1}\t{2}\t.\t+\t.\tgene_id "nameless_gene"; transcript_id "nameless_gene"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "nameless_gene"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "nameless_gene"; level 2; havana_gene "nameless_gene";'''.format(locus.chrom, locus.pos - 10000, locus.pos + 10000))

    # Transform genotype data to the right format
    
    
    # Transform genotype annotations to the right format
    

    ###########################################################
    # Then run MetaXcan using the weights database and the GWAS sumstats
    ###########################################################

    gwas_base = "/".join(locus.gwas_file.split("/")[:-1])
    gwas_file = locus.gwas_file.split("/")[-1]

    command = '''/users/mgloud/software/Metaxcan/MetaXcan/software/MetaXcan.py \
            --model_db_path {4} \
            --covariance {3} \
            --gwas_folder {2} \
            --gwas_file_pattern "{1}" \
            --snp_column SNP \
            --effect_allele_column alt \
            --non_effect_allele_column ref \
            --beta_column beta \
            --pvalue_column pvalue \
            --output_file {0}/metaxcan/metaxcan_results.txt'''.format(locus.tmpdir, gwas_file, gwas_base, covariance_file, database_file)

    ######################
    
    with open("{0}/metaxcan/metaxcan_results.txt".format(locus.tmpdir)) as f:
        f.readline()
        data = f.readline().strip().split(",")
        metaxcan_p = float(data[4])
        snps_tested = int(data[9])

    # Add results to the desired file
    with open("{0}/{1}_metaxcan_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, snps_tested, -1*math.log10(metaxcan_p)))
    
    return metaxcan_p

