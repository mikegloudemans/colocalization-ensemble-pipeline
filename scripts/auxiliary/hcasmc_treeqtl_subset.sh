join -1 1 -2 5 <(sort -k1,1 /srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_merged/treeQTL/eGenes.tsv) <(sort -k5,5 /srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gencode.v19.genes.v6p.hg19.bed) | awk '{print $8}' > /users/mgloud/projects/brain_gwas/data/gene_lists/hcasmc/treeqtl_subset.txt
