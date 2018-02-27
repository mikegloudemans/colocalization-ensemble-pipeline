#!/bin/bash

# sQTL data
sqtl_dir=/users/mgloud/projects/brain_gwas/data/eqtls/rpe/sqtls
header="feature\tchr\tsnp_pos\tref\talt\tbuild\tfeature_distance\tpvalue\tbeta\tse"
echo -e $header > $sqtl_dir/rpe_glucose_sqtls.txt
echo -e $header > $sqtl_dir/rpe_galactose_sqtls.txt

zcat /srv/persistent/bliu2/rpe/processed_data/sqtl/fastQTL/nominal/glucose/all.nominal.txt.gz | sed s/chr//g | sed s/_/\\t/g | grep -v nan | awk '{print $1 "_" $2 "\t" $0}' | cut -f2,3 --complement | sort -k2,2 -k3,3n >> $sqtl_dir/rpe_glucose_sqtls.txt
bgzip -f $sqtl_dir/rpe_glucose_sqtls.txt
tabix -S 1 -s 2 -b 3 -e 3 $sqtl_dir/rpe_glucose_sqtls.txt.gz

zcat /srv/persistent/bliu2/rpe/processed_data/sqtl/fastQTL/nominal/galactose/all.nominal.txt.gz | sed s/chr//g | sed s/_/\\t/g | grep -v nan | awk '{print $1 "_" $2 "\t" $0}' | cut -f2,3 --complement | sort -k2,2 -k3,3n >> $sqtl_dir/rpe_galactose_sqtls.txt
bgzip -f $sqtl_dir/rpe_sqtls.txt
tabix -S 1 -s 2 -b 3 -e 3 $sqtl_dir/rpe_galactose_sqtls.txt.gz


