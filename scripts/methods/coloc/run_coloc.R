suppressWarnings(suppressMessages(require(coloc)))
suppressWarnings(suppressMessages(require(readr)))

# Get input and output locations as command-line arguments
args <- commandArgs(TRUE)
infile = args[1]
N_source = as.numeric(args[2])
s_source = as.numeric(args[3])
type_source = args[4]
N_lookup = as.numeric(args[5])
s_lookup = as.numeric(args[6])
type_lookup = args[7]

# Run coloc
#data = read.table(infile, header=TRUE, sep='\t', fill=TRUE)
data = suppressWarnings(suppressMessages(as.data.frame(read_delim(infile, delim='\t'))))

# remove incomplete cases
data = data[complete.cases(data),]

# GWAS data
dataset1 = list(pvalues=data$pvalue_source, MAF=data$ref_af, s=s_source, N=N_source, type=type_source)
# eQTL data
dataset2 = list(pvalues=data$pvalue_lookup, MAF=data$ref_af, N=N_lookup, type=type_lookup)

sink(file="/dev/null")
results = coloc.abf(dataset1, dataset2)
sink()

#summary is a vector giving the number of SNPs analysed, 
#and the posterior probabilities of H0 (no causal variant), 
#H1 (causal variant for trait 1 only), 
#H2 (causal variant for trait 2 only), 
#H3 (two distinct causal variants) and 
#H4 (one common causal variant)

chr = data$chr_source[1]
feature = data$feature[1]
pos = data$seed_pos[1]
pvalue_source = min(data$pvalue_source)
pvalue_lookup = min(data$pvalue_lookup)
file_source = strsplit(strsplit(infile, "TRAIT1-")[[1]][2], "\\.TRAIT2")[[1]][1]
file_lookup = strsplit(strsplit(infile, "TRAIT2-")[[1]][2], "\\/locus")[[1]][1]

h0 = results$summary[2]
h1 = results$summary[3]
h2 = results$summary[4]
h3 = results$summary[5]
h4 = results$summary[6]

cat(chr, pos, file_source, file_lookup, pvalue_source, pvalue_lookup, feature, h0, h1, h2, h3, h4, sep='\t')

