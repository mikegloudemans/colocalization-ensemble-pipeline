#!/bin/bash

# Get GTEx eQTL files and index with tabix

cp /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_data/MatrixEQTL/allCisSNPGenePairs/* /mnt/lab_data/montgomery/mgloud/brain_gwas/data/eqtls/tabix

for f in `ls /users/mgloud/projects/brain_gwas/data/eqtls/tabix/*cis.eqtl.gz`
do
	gunzip $f
	base=${f::-3}
	cat <(echo "chr	snp_pos	ref	alt	genome	gene	beta	t-stat	pvalue") <(tail -n +2 $base | sed 's/_/\t/g' | sort -k1,1n -k2,2n) | bgzip > $base.gz
	tabix $base.gz -b 2 -e 2 -s 1 -S 1 
	rm $base
done

