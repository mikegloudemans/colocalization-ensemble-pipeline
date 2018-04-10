#!/bin/bash
cat <(zcat /mnt/lab_data/montgomery/shared/1KG/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | tail -n +253 | head -n 1) <(zcat /mnt/lab_data/montgomery/shared/1KG/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk '{if (($3 == "rs3138142") || ($3 == "rs3138141")) print $0}') | column -t > 1kg_alleles.txt
