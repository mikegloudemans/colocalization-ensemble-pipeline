zcat ../../data/eqtls/hcasmc/hcasmc_eqtls.txt.gz | tail -n +2 | awk '{if ($26 < 0.00005) print $1}' | sort | uniq > ../../data/gene_lists/hcasmc/eqtl_subset_5e-5.txt
