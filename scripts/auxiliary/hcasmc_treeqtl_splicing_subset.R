require(dplyr)
require(readr)
require(qvalue)

d = read_delim(gzfile("/srv/persistent/bliu2/HCASMC_eQTL/processed_data/sqtl/fastQTL/permutation/all.permutation.txt.gz"), delim="\t", col_names=FALSE)
d$qvalue = qvalue(d$X16)$qvalue

print(dim(d))
d = d[d$qvalue < 0.05,]
print(dim(d))

write_tsv(as.data.frame(d[,1]), "/users/mgloud/projects/brain_gwas/data/gene_lists/hcasmc_treeqtl_sqtls.txt", col_names=FALSE)
