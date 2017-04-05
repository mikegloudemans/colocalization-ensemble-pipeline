# Author: Mike Gloudemans
# Revised 4/5/2017

# Analyze distribution of candidate colocalization events after
# colocalization detection, to assess the extent of sensitivity
# to tissue sample sizes.

library(dplyr)

coloc_threshold = 0.01

# Load colocalization events table, and tissue sample sizes.
d = read.table("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/scz2_snp_results_txt_finemap_clpp_status.txt", header = FALSE)
dTissues = read.csv("/mnt/lab_data/montgomery/joed3/GTExCisEqtls/data/gtex.v6p.egenes.summary.txt", sep="\t")
colnames(d) <- c("snp", "Tissue", "gene", "conditional_level", "num_tested", "clpp")
d$coloc <- d$clpp > coloc_threshold

tissue_coloc = d %>% group_by(Tissue) %>% summarize(num_colocs=sum(coloc), pct_coloc = sum(coloc) / length(coloc) * 100)


# Plot distribution of number of colocalization events per tissue.
png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/finemap_colocs_per_tissue.png", units="in", width=5, height=5, res=360)
hist(tissue_coloc$num_colocs)
dev.off()

# Plot relationship between number of colocalization events and
# sample size.
merged = merge(dTissues, tissue_coloc, by="Tissue") 
brain_merged = merged[grep("Brain", merged$Tissue),] 

png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/finemap_colocs_vs_sample_size.png", units="in", width=5, height=5, res=360)

	plot(merged$cisN, merged$num_colocs, xlab="Sample size", ylab="Number of colocs (CLPP > 0.01)", pch=19)
	points(brain_merged$cisN, brain_merged$num_colocs, pch=19, col="gold")
dev.off()

# Get total number of colocalization events at each gene with
# at least one colocalization. Then plot the distribution.
gene_coloc_counts = d %>% group_by(gene) %>% summarize(num_colocs = sum(coloc))
# Fraction of tested genes that had at least one coloc event
sum(gene_coloc_counts$num_colocs >= 1) / length(gene_coloc_counts$num_colocs)
gene_coloc_counts = gene_coloc_counts[gene_coloc_counts$num_colocs >= 1,]

png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/colocs_per_gene.png", units="in", width=5, height=5, res=360)
	hist(gene_coloc_counts$num_colocs, xlab="Number of colocalization events", ylab="Number of genes")
dev.off()

# Get total number of colocalization events at each gene with
# at least one colocalization. Then plot the distribution.
snp_coloc_counts = d %>% group_by(snp) %>% summarize(num_colocs = sum(coloc))
# Fraction of tested SNPs that had at least one coloc event
# NOTE: This number could be a bit inflated since some SNPs
# were likely not tested at all if they had no nearby genes
sum(snp_coloc_counts$num_colocs >= 1) / length(snp_coloc_counts$num_colocs)
snp_coloc_counts = snp_coloc_counts[snp_coloc_counts$num_colocs >= 1,]

png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/colocs_per_snp.png", units="in", width=5, height=5, res=360)
	hist(snp_coloc_counts$num_colocs, xlab="Number of colocalization events", ylab="Number of snps")
dev.off()

# Get total number of colocalization events at each gene with
# at least one colocalization. Then plot the distribution.
paired_coloc_counts = d %>% group_by(snp, gene) %>% summarize(num_colocs = sum(coloc))
# Fraction of tested SNP-gene pairs that had at least one coloc event
sum(paired_coloc_counts$num_colocs >= 1) / length(paired_coloc_counts$num_colocs)
paired_coloc_counts = paired_coloc_counts[paired_coloc_counts$num_colocs >= 1,]

png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/colocs_per_snp_gene_pair.png", units="in", width=5, height=5, res=360)
        hist(paired_coloc_counts$num_colocs, xlab="Number of colocalization events", ylab="Number of snp-gene pairs")
dev.off()

# Split data into those colocalization events shared across tissues, and those
# unique to a single tissue.
shared_coloc = filter(paired_coloc_counts, num_colocs > 1) 
unique_coloc = filter(paired_coloc_counts, num_colocs == 1) 
dim(unique_coloc)
dim(shared_coloc)

# Plot relationship between number of unique colocalization events 
# and tissue sample size.
colocs_only = filter(d, coloc==TRUE)
unique_colocs_only = merge(unique_coloc, colocs_only, by=c("snp", "gene"))
unique_coloc_counts = unique_colocs_only %>% group_by(Tissue) %>% summarize(num_colocs = length(coloc))
unique_coloc_counts = merge(unique_coloc_counts, dTissues, by="Tissue")
brain_unique_counts = unique_coloc_counts[grep("Brain", unique_coloc_counts$Tissue),]

png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/finemap_unique_colocs_vs_sample_size.png", units="in", width=5, height=5, res=360)
	plot(unique_coloc_counts$cisN, unique_coloc_counts$num_colocs, xlab="Sample size", ylab="Number of unique colocalization events", pch=19)
	points(brain_unique_counts$cisN, brain_unique_counts$num_colocs, pch=19, col="gold")
dev.off()

cor.test(unique_coloc_counts$cisN, unique_coloc_counts$num_colocs)

# Plot relationship between number of shared colocalization events 
# and tissue sample size.
shared_colocs_only = merge(shared_coloc, colocs_only, by=c("snp", "gene"))
shared_coloc_counts = shared_colocs_only %>% group_by(Tissue) %>% summarize(num_colocs = length(coloc))
shared_coloc_counts = merge(shared_coloc_counts, dTissues, by="Tissue")
brain_shared_counts = shared_coloc_counts[grep("Brain", shared_coloc_counts$Tissue),]

png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/finemap_shared_colocs_vs_sample_size.png", units="in", width=5, height=5, res=360)
	plot(shared_coloc_counts$cisN, shared_coloc_counts$num_colocs, xlab="Sample size", ylab="Number of shared colocalization events", pch=19)
	points(brain_shared_counts$cisN, brain_shared_counts$num_colocs, pch=19, col="gold")
dev.off()

cor.test(shared_coloc_counts$cisN, shared_coloc_counts$num_colocs)

