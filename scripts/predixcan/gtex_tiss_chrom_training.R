# modified to read in file paths as arguments

# setwd("/group/im-lab/nas40t2/scott/gtex_v7_imputed_europeans/model_training/scripts/") # modified, commented
source("gtex_v7_nested_cv_elnet.R")
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
tiss <- argv[1]
chrom <- argv[2]

# snp_annot_file <- "../../prepare_data/genotype/gtex_v7_hapmapceu_snp_annot.chr" %&% chrom %&% ".txt"
snp_annot_file <- argv[3]
# gene_annot_file <- "../../prepare_data/expression/gencode.v19.genes.patched_contigs.parsed.txt"
gene_annot_file <- argv[4]
# genotype_file <- "../../prepare_data/genotype/gtex_v7_eur_imputed_maf_0.01_R2_0.8_chr" %&% chrom %&% ".txt.gz"
genotype_file <- argv[5]
# expression_file <- "../../prepare_data/expression/" %&% tiss %&% "_Analysis.expression.txt"
expression_file <- argv[6]
# covariates_file <- "../../prepare_data/covariates/" %&% tiss %&% "_Analysis.combined_covariates.txt"
covariates_file <- argv[7]
prefix <- tiss %&% "_nested_cv"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)


