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
output_dir = args[2]

if (length(args) > 2)
{
	eqtl_filter = args[3]
}

# Load data
results = read_delim(input_file, col_names=FALSE, delim ="\t")
colnames(results) = c("snp", "eqtl_study", "gwas_study", "gene", "conditional_level", "snps_tested", "clpp_score", "gwas_log_pval", "eqtl_log_pval")
chrs = read_delim("/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.dict", delim="\t", skip=1, col_names=FALSE)

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

# Temporary changes to investigate alternative plot scales 
#results$clpp_score = sapply(results$clpp_score, function(x){min(0.02, x)})

# Getting chromosomal location for each data point
full_results = inner_join(results, chrs, by=c("chrom"))
full_results$abs_pos = full_results$cumsum + full_results$pos

all_endpoints = c(chrs$cumsum, tail(chrs$cumsum,1)+tail(chrs$length,1))
midpoints = sapply(1:(length(all_endpoints)-1), function(x){(all_endpoints[x] + all_endpoints[x+1]) / 2})

rich6equal = c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300")
palette(rich6equal)

for (gwas_file in unique(full_results$gwas_study))
{
	for (eqtl_file in unique(full_results$eqtl_study))
	{
		results_subset = full_results[(full_results$gwas_study == gwas_file) & (full_results$eqtl_study == eqtl_file),]

		top_results = results_subset[results_subset$clpp_score > 0.005,]
		if (dim(top_results)[1] == 0)
		{
			next
		}
		# color scheme from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
		png(paste(output_dir, "/", gwas_file, "_", eqtl_file, ".png", sep=""), units="in",width=15, height=4, res=360)
			plot(results_subset$abs_pos, results_subset$clpp_score,pch=19, xlab="Chromosome", ylab="CLPP score", col=results_subset$chrom, cex=1.5, xaxt='n')
			with(top_results, text(clpp_score~abs_pos, labels = gene, pos = 2))
			axis(side = 1, at = midpoints, tick=FALSE, labels=1:22)
		dev.off()
		write_delim(results_subset, paste(output_dir, "/", gwas_file, "_", eqtl_file, ".txt", sep=""), delim="\t", col_names=TRUE)
	}
}

