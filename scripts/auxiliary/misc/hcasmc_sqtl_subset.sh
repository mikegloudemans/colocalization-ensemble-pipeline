zcat ../../data/eqtls/hcasmc/sqtls/hcasmc_sqtls.txt.gz | tail -n +2 | awk '{if ($9 < 0.00005) print $6}' | sort | uniq > ../../data/gene_lists/hcasmc/sqtl_subset_5e-5.txt
