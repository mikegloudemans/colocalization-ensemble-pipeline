#!/bin/bash
# Author: Mike Gloudemans
# 9/21/2017
# Run BEAGLE to phase DGN genotypes

for i in {22..22}; do


	zcat /mnt/lab_data/montgomery/shared/datasets/dgn/genotypes_from_kim/imputation/otherFormats/filtered/filtered-v2/GenRED.II.autosomal.final.chr$i.QCTOOL.v2.fil-v2.recode.sorted.vcf.gz | awk '{if(($4 != "-") && ($5 != "-")) print $0}' > dgn_tmp$i.vcf
	python filter_maf.py /mnt/lab_data/montgomery/shared/1KG/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 1kgenomes_tmp$i.vcf

	# Need to determine whether missing genotypes should be imputed or not
	java -Xss5m -Xmx8g -jar /users/mgloud/software/beagle/beagle.08Jun17.d8b.jar gt=dgn_tmp$i.vcf out=chr$i ref=1kgenomes_tmp$i.vcf impute=true nthreads=8 ibd=false window=10000 overlap=2000

	rm dgn_tmp$i.vcf
	rm 1kgenomes_tmp$i.vcf
done
