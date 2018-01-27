#!/bin/bash

# Date created: 11/10/2017
# Author: Mike Gloudemans

for tissue in `ls /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/all_associations/`; do
	tissue=`echo $tissue | sed s/\.gz//g`
	header="gene\tchr\tsnp_pos\tref\talt\tbuild\ttss_distance\tma_samples\tma_count\tmaf\tpvalue\tbeta\tse"
	echo -e $header > /users/mgloud/projects/brain_gwas/data/eqtls/gtex_v7/$tissue
	zcat /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/all_associations/$tissue.gz | tail -n +2 | sed s/_/\\t/g | grep -v nan | sort -k2,2 -k3,3n >> /users/mgloud/projects/brain_gwas/data/eqtls/gtex_v7/$tissue
	bgzip -f /users/mgloud/projects/brain_gwas/data/eqtls/gtex_v7/$tissue
	tabix -f -S 1 -s 2 -b 3 -e 3 /users/mgloud/projects/brain_gwas/data/eqtls/gtex_v7/$tissue.gz
done
