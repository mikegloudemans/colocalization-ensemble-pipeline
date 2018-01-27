#!/bin/bash

# Note: I think ref/alt are oriented correctly. See
# the 23andme documentation for all GWAS analyses to confirm.
# It depends which one we consider the "B" allele; I'd assume
# it's whichever one is listed second.
echo -e "rsid\tchr\tsnp_pos\tref\talt\tpvalue\tbeta\tse\tpass" > /users/mgloud/projects/brain_gwas/data/gwas/rpe/23andme_myopia.prepared.txt
join -t $'\t' <(tail -n +2 /users/mgloud/projects/brain_gwas/data/gwas/rpe/23andme/4.1_annotation/gwas_anno-4.1/all_snp_info-4.1.txt | cut -f1,4,5,6,7 | sed s/\\//\\t/g) <(tail -n +2 /users/mgloud/projects/brain_gwas/data/gwas/rpe/23andme/nearsightedness-4.1.dat | cut -f1,3,4,5,6) | cut -f1 --complement | awk '{if ($9 == "Y") print $0}' | sort -k2,2 -k3,3n >> /users/mgloud/projects/brain_gwas/data/gwas/rpe/23andme_myopia.prepared.txt

bgzip -f /users/mgloud/projects/brain_gwas/data/gwas/rpe/23andme_myopia.prepared.txt
tabix -f -s 2 -b 3 -e 3 -S 1 /users/mgloud/projects/brain_gwas/data/gwas/rpe/23andme_myopia.prepared.txt.gz
