require(coloc)

# Get input and output locations as command-line arguments
args <- commandArgs(TRUE)
infile = args[1]
N_gwas = as.numeric(args[2])
s_gwas = as.numeric(args[3])
type_gwas = args[4]
N_eqtl = as.numeric(args[5])

# Run coloc
data = read.table(infile, header=TRUE)

# GWAS data
#dataset1 = list(pvalues=data$pvalue_gwas, MAF=data$allele_freq, s=s_gwas, N=N_gwas, type=type_gwas)
dataset1 = list(pvalues=data$pvalue_gwas, MAF=data$Kgenomes_maf, s=s_gwas, N=N_gwas, type=type_gwas)
# eQTL data
#dataset2 = list(pvalues=data$pvalue_eqtl, MAF=data$allele_freq, N=N_eqtl, type="quant")
dataset2 = list(pvalues=data$pvalue_eqtl, MAF=data$Kgenomes_maf, N=N_eqtl, type="quant")

sink(file="/dev/null")
results = coloc.abf(dataset1, dataset2)
sink()

cat(results$summary[6])
