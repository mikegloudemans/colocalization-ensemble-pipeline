suppressMessages(require(coloc))

# Get input and output locations as command-line arguments
args <- commandArgs(TRUE)
infile = args[1]
N_gwas = as.numeric(args[2])
s_gwas = as.numeric(args[3])
type_gwas = args[4]
N_eqtl = as.numeric(args[5])

# Run coloc
data = read.table(infile, header=TRUE, sep=',', fill=TRUE)

# remove incomplete cases
data = data[complete.cases(data),]

# GWAS data
#dataset1 = list(pvalues=data$pvalue_gwas, MAF=data$allele_freq, s=s_gwas, N=N_gwas, type=type_gwas)
dataset1 = list(pvalues=data$pvalue_gwas, MAF=data$effect_af_gwas, s=s_gwas, N=N_gwas, type=type_gwas)
# eQTL data
#dataset2 = list(pvalues=data$pvalue_eqtl, MAF=data$allele_freq, N=N_eqtl, type="quant")
dataset2 = list(pvalues=data$pvalue_eqtl, MAF=data$effect_af_eqtl, N=N_eqtl, type="quant")

sink(file="/dev/null")
results = coloc.abf(dataset1, dataset2)
sink()

#summary is a vector giving the number of SNPs analysed, 
#and the posterior probabilities of H0 (no causal variant), 
#H1 (causal variant for trait 1 only), 
#H2 (causal variant for trait 2 only), 
#H3 (two distinct causal variants) and 
#H4 (one common causal variant)

h0 = results$summary[2]
h1 = results$summary[3]
h2 = results$summary[4]
h3 = results$summary[5]
h4 = results$summary[6]

cat(h0, h1, h2, h3, h4, sep='\t')

