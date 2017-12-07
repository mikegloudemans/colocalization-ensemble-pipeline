#!/bin/bash

# Author: Mike Gloudemans
#
# Commands to wrangle the RPE data into a usable format for coloc running

# GWAS data

header="Marker\tchr\tsnp_pos\tref\talt\tNcases\tNcontrols\tpvalue\tdirection"

# Note that we need to swap columns 4 and 5 because the effect allele is listed first in the original data.
echo -e $header > /users/mgloud/projects/brain_gwas/data/gwas/rpe/Fritsche_sorted.txt
tail -n +2 /srv/persistent/bliu2/rpe/data/gwas/Fritsche_2015_AdvancedAMD.txt | sort -k2,2 -k3,3n |  \
	awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $4 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' >> /users/mgloud/projects/brain_gwas/data/gwas/rpe/Fritsche_sorted.txt
bgzip -f /users/mgloud/projects/brain_gwas/data/gwas/rpe/Fritsche_sorted.txt
tabix -S 1 -s 2 -b 3 -e 3 /users/mgloud/projects/brain_gwas/data/gwas/rpe/Fritsche_sorted.txt.gz

# eQTL data
header="gene\trsid\tchr\tsnp_pos\tref\talt\tallele_freq\thwe_chisq\tIA\tqval\tchisq\tpi\tmapping_error_rate\tref_bias\toverdispersion\tsnp_regional_id\tfeature_snp_count\ttested_snp_count\tnull_iterations\talt_iterations\trandom_tie_location\tlogl_null\tconvergence\tr2_fSNP\tr2_rSNP"

cat <(echo -e $header) <(cat /users/mgloud/projects/brain_gwas/data/eqtls/output/glucose/joint/*/* | sed 's/chr//g') | sort -k3,3n -k4,4n | bgzip > /users/mgloud/projects/brain_gwas/data/eqtls/output/glucose.eqtls.txt.gz
tabix -p bed -S 1 -s 3 -b 4 /users/mgloud/projects/brain_gwas/data/eqtls/output/glucose.eqtls.txt.gz

cat <(echo -e $header) <(cat /users/mgloud/projects/brain_gwas/data/eqtls/output/galactose/joint/*/* | sed 's/chr//g') | sort -k3,3n -k4,4n | bgzip > /users/mgloud/projects/brain_gwas/data/eqtls/output/galactose.eqtls.txt.gz
tabix -p bed -S 1 -s 3 -b 4 /users/mgloud/projects/brain_gwas/data/eqtls/output/galactose.eqtls.txt.gz

# sQTL data
sqtl_dir=/users/mgloud/projects/brain_gwas/data/eqtls/rpe/sqtls
header="feature\tchr\tsnp_pos\tref\talt\tbuild\tfeature_distance\tpvalue\tbeta\tse"
echo -e $header > $sqtl_dir/rpe_sqtls.txt

zcat /srv/persistent/bliu2/rpe/processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz | sed s/chr//g | sed s/_/\\t/g | grep -v nan | awk '{print $1 "_" $2 "\t" $0}' | cut -f2,3 --complement | sort -k2,2 -k3,3n >> $sqtl_dir/rpe_sqtls.txt

bgzip -f $sqtl_dir/rpe_sqtls.txt
tabix -S 1 -s 2 -b 3 -e 3 $sqtl_dir/rpe_sqtls.txt.gz

