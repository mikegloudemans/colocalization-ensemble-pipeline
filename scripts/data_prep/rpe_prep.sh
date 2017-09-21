#!/bin/bash

# Author: Mike Gloudemans
#
# Commands to wrangle the RPE data into a usable format for coloc running

#cp -r /srv/persistent/bliu2/rpe/processed_data/rasqual/output/ /users/mgloud/projects/brain_gwas/data/eqtls/

header="gene\trsid\tchr\tsnp_pos\tref\talt\tallele_freq\thwe_chisq\tIA\tqval\tchisq\teffect_size\tmapping_error_rate\tref_bias\toverdispersion\tsnp_regional_id\tfeature_snp_count\ttested_snp_count\tnull_iterations\talt_iterations\trandom_tie_location\tlogl_null\tconvergence\tr2_fSNP\tr2_rSNP"

cat <(echo -e $header) <(cat /users/mgloud/projects/brain_gwas/data/eqtls/output/glucose/joint/*/* | sed 's/chr//g') | sort -k3,3n -k4,4n | bgzip > /users/mgloud/projects/brain_gwas/data/eqtls/output/glucose.eqtls.txt.gz
tabix -p bed -S 1 -s 3 -b 4 /users/mgloud/projects/brain_gwas/data/eqtls/output/glucose.eqtls.txt.gz

cat <(echo -e $header) <(cat /users/mgloud/projects/brain_gwas/data/eqtls/output/galactose/joint/*/* | sed 's/chr//g') | sort -k3,3n -k4,4n | bgzip > /users/mgloud/projects/brain_gwas/data/eqtls/output/galactose.eqtls.txt.gz
tabix -p bed -S 1 -s 3 -b 4 /users/mgloud/projects/brain_gwas/data/eqtls/output/galactose.eqtls.txt.gz
