#!/usr/bin/Rscript
# Author: Mike Gloudemans
# Date created: 12/7/2017
# Plot genome-wide colocalization scores for an analysis.

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

args <- commandArgs(TRUE)
input_file = args[1]
output_file = args[2]
eqtl_subset= ""
if (!is.na(args[3]))
{
	eqtl_subset = args[3]
}

# Load data
#results = read_delim("../output/2017-11-13_08-51-41_hcasmc_full/UKBB_GWAS1KG_EXOME_CAD_SOFT_META_PublicRelease_300517_txt_gz_finemap_clpp_status.txt", col_names=FALSE, delim ="\t")
results = read_delim(input_file, col_names=FALSE, delim ="\t")
colnames(results) = c("snp", "eqtl_study", "gwas_study", "gene", "conditional_level", "snps_tested", "clpp_score")
chrs = read_delim("/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.dict", delim="\t", skip=1, col_names=FALSE)

if (eqtl_subset != "")
{
	results = results[results$eqtl_study == eqtl_subset,]
}

# Munging to get cumulative position in genome for each chromosome
chrs$chrom = as.numeric(sapply(chrs[,2], function(x){substring(x,7)}))
chrs = chrs[!is.na(chrs$chrom),]
chrs$length = as.numeric(sapply(chrs[,3], function(x){as.numeric(substring(x,4))}))
chrs = chrs[order(chrs$chrom),]
chrs$cumsum = sapply(1:dim(chrs)[1], function(x){sum(chrs$length[(1:x)-1])})	# Note unorthodox parentheses placement in indexing vector; this is essential
chrs = chrs[,-c(1:5)]

# Munging results data
results$chrom = as.numeric(sapply(unlist(results[,1]), function(x){strsplit(x,"_")[[1]][1]}))
results$pos = as.numeric(sapply(unlist(results[,1]), function(x){strsplit(x,"_")[[1]][2]}))

# Getting chromosomal location for each data point
full_results = inner_join(results, chrs, by=c("chrom"))
full_results$abs_pos = full_results$cumsum + full_results$pos

all_endpoints = c(chrs$cumsum, tail(chrs$cumsum,1)+tail(chrs$length,1))
midpoints = sapply(1:(length(all_endpoints)-1), function(x){(all_endpoints[x] + all_endpoints[x+1]) / 2})

# color scheme from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
png(output_file, units="in",width=15, height=4, res=360)
	rich6equal = c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300")
	palette(rich6equal)
	plot(full_results$abs_pos, full_results$clpp_score,pch=19, xlab="Chromosome", ylab="CLPP score", col=full_results$chrom, cex=1.5, xaxt='n')
	top_results = full_results[full_results$clpp_score > 0.005,]
	with(top_results, text(clpp_score~abs_pos, labels = gene, pos = 2))
	axis(side = 1, at = midpoints, tick=FALSE, labels=1:22)
dev.off()
