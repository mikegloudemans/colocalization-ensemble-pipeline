#join -1 1 -2 3 <(sort rsid_list.txt) <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) > snps_with_hg19_positions.txt
join -1 1 -2 3 <(sort insulin_rsids.txt) <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) | grep -v hap > insulin_rsids_with_hg19_positions.txt

