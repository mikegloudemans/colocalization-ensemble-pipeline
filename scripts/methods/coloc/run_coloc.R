suppressWarnings(suppressMessages(require(coloc)))
suppressWarnings(suppressMessages(require(readr)))

# Get input and output locations as command-line arguments
args <- commandArgs(TRUE)
infile = args[1]
N_trait1 = as.numeric(args[2])
s_trait1 = as.numeric(args[3])
type_trait1 = args[4]
N_trait2 = as.numeric(args[5])
s_trait2 = as.numeric(args[6])
type_trait2 = args[7]

# Run coloc
#data = read.table(infile, header=TRUE, sep='\t', fill=TRUE)
data = suppressWarnings(suppressMessages(as.data.frame(read_delim(infile, delim='\t'))))

# remove incomplete cases
#data = data[complete.cases(data),]

# GWAS data
dataset1 = list(pvalues=data$pvalue_trait1, MAF=data[["ref_af_vcf1"]], s=s_trait1, N=N_trait1, type=type_trait1)
# eQTL data
dataset2 = list(pvalues=data$pvalue_trait2, MAF=data[["ref_af_vcf1"]], N=N_trait2, type=type_trait2)

sink(file="/dev/null")
results = coloc.abf(dataset1, dataset2)
sink()

#summary is a vector giving the number of SNPs analysed, 
#and the posterior probabilities of H0 (no causal variant), 
#H1 (causal variant for trait 1 only), 
#H2 (causal variant for trait 2 only), 
#H3 (two distinct causal variants) and 
#H4 (one common causal variant)

N_snps = dim(data)[1]

chr = data$chr_trait1[1]
feature1 = data$feature_trait1[1]
feature2 = data$feature_trait2[1]
pos = data$seed_pos[1]
pvalue_trait1 = min(data$pvalue_trait1)
pvalue_trait2 = min(data$pvalue_trait2)
file_trait1 = strsplit(strsplit(infile, "TRAIT1-")[[1]][2], "\\.TRAIT2")[[1]][1]
file_trait2 = strsplit(strsplit(infile, "TRAIT2-")[[1]][2], "\\/locus")[[1]][1]

h0 = results$summary[2]
h1 = results$summary[3]
h2 = results$summary[4]
h3 = results$summary[5]
h4 = results$summary[6]

cat(chr, pos, file_trait1, file_trait2, pvalue_trait1, pvalue_trait2, feature1, feature2, N_snps, h0, h1, h2, h3, h4, sep='\t')
