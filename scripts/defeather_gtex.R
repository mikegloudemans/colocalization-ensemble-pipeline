# Switch GTEx gene expression data from a feather format to a TSV format,
# which is much easier to deal with in Python currently

require(feather)
require(readr)

data = read_feather("data/gtex/GTEx_Analysis_2017-06-05_v8/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.ft")
write_delim(data, "data/gtex-reformatted/eqtl.tpm.txt", delim="\t", col_names=TRUE)
