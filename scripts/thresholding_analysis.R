# Analyze distribution of candidate colocalization events after thresholding
# method is applied, to get a sense of whether this method is less biased
# by tissue sample sizes.

# Load colocalization events table, and tissue sample sizes.
d = read.table("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/scz2_snp_results_txt_finemap_coloc_status.txt", header = TRUE)
dTissues = read.csv("/users/joed3/GTExCisEqtls/data/gtex.v6p.egenes.summary.txt", sep="\t")

# Plot distribution of number of colocalization events per tissue.
png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-06_by_1.0e-07/coloc_status/plots/colocs_per_tissue.png", units="in", width=5, height=5, res=360)
	hist(colSums(d[,3:dim(d)[2]]), xlab="Number of colocalization events")
dev.off()

# Plot relationship between number of colocalization events and
# sample size. 
coloc_counts <- colSums(d[,3:dim(d)[2]])
coloc_counts <- as.data.frame(cbind(names(coloc_counts),coloc_counts))
coloc_counts$coloc_counts <- as.numeric(as.character(coloc_counts$coloc_counts))
names(coloc_counts)[1] <- "Tissue"
tissues_with_counts = merge(dTissues, coloc_counts, by="Tissue")

brain_counts = tissues_with_counts[grep("Brain", tissues_with_counts$Tissue),]

png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-06_by_1.0e-07/coloc_status/plots/colocs_vs_sample_size.png", units="in", width=5, height=5, res=360)
	plot(tissues_with_counts$Size, tissues_with_counts$coloc_counts, xlab="Sample size", ylab="Number of candidate colocalization events", pch=19)
	points(brain_counts$Size, brain_counts$coloc_counts, pch=19, col="gold")
dev.off()

# Get total number of tissues colocalizing at each gene/SNP pair with
# at least one colocalization. Then plot the distribution.
event_count <- rowSums(d[,3:dim(d)[2]])
event_count <- event_count[event_count != 0]
png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-06_by_1.0e-07/coloc_status/plots/colocs_per_snp_gene_pair.png", units="in", width=5, height=5, res=360)
	hist(event_count, xlab="Number of colocalization events")
dev.off()

# Subset colocalization table down to the rows that have at least
# one candidate colocalization event.
coloc_events = d[rowSums(d[,3:dim(d)[2]]) != 0,]

# Split data into those colocalization events shared across tissues, and those
# unique to a single tissue.
shared_coloc <- coloc_events[rowSums(coloc_events[,3:dim(coloc_events)[2]]) > 15,]
unique_coloc <- coloc_events[rowSums(coloc_events[,3:dim(coloc_events)[2]]) <= 15,]
dim(unique_coloc)
dim(shared_coloc)

# Plot relationship between number of unique colocalization events 
# and tissue sample size.
unique_coloc_counts <- colSums(unique_coloc[,3:dim(unique_coloc)[2]])
unique_coloc_counts <- as.data.frame(cbind(names(unique_coloc_counts),unique_coloc_counts))
unique_coloc_counts$coloc_counts <- as.numeric(as.character(unique_coloc_counts$unique_coloc_counts))
names(unique_coloc_counts)[1] <- "Tissue"
tissues_with_unique_counts = merge(dTissues, unique_coloc_counts, by="Tissue")

brain_unique_counts = tissues_with_unique_counts[grep("Brain", tissues_with_unique_counts$Tissue),]

png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-06_by_1.0e-07/coloc_status/plots/unique_colocs_vs_sample_size.png", units="in", width=5, height=5, res=360)
	plot(tissues_with_unique_counts$Size, tissues_with_unique_counts$unique_coloc_counts, xlab="Sample size", ylab="Number of candidate colocalization events", pch=19)
	points(brain_unique_counts$Size, brain_unique_counts$unique_coloc_counts, pch=19, col="gold")
dev.off()

# Plot relationship between number of shared colocalization events 
# and tissue sample size.
shared_coloc_counts <- colSums(shared_coloc[,3:dim(shared_coloc)[2]])
shared_coloc_counts <- as.data.frame(cbind(names(shared_coloc_counts),shared_coloc_counts))
shared_coloc_counts$shared_coloc_counts <- as.numeric(as.character(shared_coloc_counts$shared_coloc_counts))
names(shared_coloc_counts)[1] <- "Tissue"
tissues_with_shared_counts = merge(dTissues, shared_coloc_counts, by="Tissue")

brain_shared_counts = tissues_with_shared_counts[grep("Brain", tissues_with_shared_counts$Tissue),]

png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-06_by_1.0e-07/coloc_status/plots/shared_colocs_vs_sample_size.png", units="in", width=5, height=5, res=360)
	plot(tissues_with_shared_counts$Size, tissues_with_shared_counts$shared_coloc_counts, xlab="Sample size", ylab="Number of candidate colocalization events", pch=19)
	points(brain_shared_counts$Size, brain_shared_counts$shared_coloc_counts, pch=19, col="gold")
dev.off()

