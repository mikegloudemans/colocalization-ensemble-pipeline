# Author: Mike Gloudemans
# Date created: 4/10/2018

require(reshape)
require(ggplot2)
require(dplyr)

#
# Load files
#

files = c("/users/mgloud/projects/brain_gwas/output/completed/eriks_grant/2018-04-19_07-49-06_insulin_for_eriks_grant/FastGlu_MAGIC_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/completed/eriks_grant/2018-04-19_07-49-06_insulin_for_eriks_grant/FastInsu_adjBMI_MAGIC_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/completed/eriks_grant/2018-04-19_07-49-06_insulin_for_eriks_grant/HDL_GCLC_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/completed/eriks_grant/2018-04-19_07-49-06_insulin_for_eriks_grant/TG_GCLC_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/completed/eriks_grant/2018-04-19_07-49-06_insulin_for_eriks_grant/WHRadjBMI_GIANT_Mixed_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/completed/eriks_grant/2018-04-23_16-16-49_insulin_for_eriks_grant_bonus_fastinsu/GWAS_Glycemic-Traits_Dupuis_2010_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/completed/eriks_grant/2018-04-23_15-04-54_insulin_for_eriks_grant_MI/GWAS_GENESIS_MI_txt_gz_finemap_clpp_status.txt"
)

# Finemap results
finemap_results = list()
for (i in 1:length(files))
{
	finemap_results[[i]] = read.table(files[i], header=TRUE)
}
finemap = do.call(rbind, finemap_results)
finemap$ensembl = sapply(as.character(finemap$feature), function(x) {strsplit(x, "\\.")[[1]][1]})
finemap$gwas_trait = gsub(".prepared.txt.gz", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("/users/mgloud/projects/gwas/data/prepared/", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("/users/mgloud/projects/gwas/data/GENESIS/formatted/", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("GWAS_", "", finemap$gwas_trait)
finemap$gwas_trait = gsub(".txt.gz", "", finemap$gwas_trait)
#finemap$gwas_trait = sapply(as.character(finemap$gwas_trait), function(x) 
#       {
#	       s = strsplit(x, "_")[[1]]
#	       print(s)
#	       return(paste(s[1:(length(s)-3)], collapse="_"))
#       })

finemap$eqtl_file = sapply(as.character(finemap$eqtl_file), function(x)
       {
               s = strsplit(x, "_")[[1]]
               return(paste(s[1:(length(s)-3)], collapse="_"))
       })

# Load rsIDs
rsids = read.table("/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/insulin_rsids_with_hg19_positions.txt", header=FALSE)
rsids[,2] = gsub("chr", "", rsids[,2])
rsids$ref_snp = paste(rsids[,2], rsids[,3], sep="_")
rsids = rsids[,-c(2,3)]
colnames(rsids)[1] = "rsid"

# Get HGNC gene names
genes = read.table("/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/mart_export.txt", header=TRUE, sep="\t")
genes = genes[,-2]
names(genes) = c("ensembl", "hgnc")

genes = genes[!duplicated(genes$ensembl),]

finemap = merge(rsids, finemap)
finemap = merge(genes, finemap)

#
# First, get the list of SNP/gene pairs passing the desired cutoff, so that we can
# show them in individual plots
#
# NOTE: Hard to determine where to put the cutoff here, because some of our GWAS have 
# very few SNPs tested, which leads to inflated CLPP scores. May be possible to correct
# that with a simple regression if needed, but let's skip that for now.
#
dim(finemap)	# Total number of tests performed
pass_pval_cutoffs = finemap[finemap$X.log_gwas_pval > 5 & finemap$X.log_eqtl_pval > 5,]
#pass_pval_cutoffs = finemap[finemap$X.log_gwas_pval > 5 & finemap$X.log_eqtl_pval > 7.3,]
dim(pass_pval_cutoffs)	# Total number of tests (gwas SNP - eqtl gene pairs) passing eQTL and GWAS pval thresholds
#pre_final_set = pass_pval_cutoffs[pass_pval_cutoffs$clpp_mod > 0.3,]
pre_final_set = pass_pval_cutoffs
dim(pre_final_set)	# Total number of tests passing our colocalization criteria
pre_final_set = pre_final_set[rev(order(pre_final_set$clpp_mod)),]

pre_final_set = pre_final_set[!duplicated(paste(pre_final_set$feature, pre_final_set$eqtl_file, pre_final_set$base_gwas_file)),]

# Remove duplicates from this list
final_set = pre_final_set[!duplicated(paste(pre_final_set$ref_snp, pre_final_set$feature, sep="_")),]

#
# Then, for each SNP-gene pair passing the test, show full summary statistics across
# all eQTL files and GWAS
#
eqtl_files = unique(finemap$eqtl_file)
gwas_traits = unique(finemap$gwas_trait)

# For each SNP-gene pair that had at least one colocalization
for (i in 1:dim(final_set)[1])
{
	# Make a matrix to show colocalization across all tissues
	my_matrix = array(0, c(length(eqtl_files), length(gwas_traits)))
	rownames(my_matrix) = eqtl_files
	colnames(my_matrix) = gwas_traits
	trunc_matrix = my_matrix	# This matrix will only include the SNPs if they pass the p-value cutoffs

	# For each eQTL file
	for (j in 1:length(eqtl_files))
	{
		# For each GWAS file
		for (k in 1:length(gwas_traits))
		{
			# Fill matrix accordingly, or fill with NA
			result = finemap[finemap$eqtl_file == eqtl_files[j] & finemap$gwas_trait == gwas_traits[k] & 
					 finemap$ref_snp == final_set$ref_snp[i] & finemap$feature == final_set$feature[i],]
			result_trunc = finemap[finemap$eqtl_file == eqtl_files[j] & finemap$gwas_trait == gwas_traits[k] & 
					 finemap$ref_snp == final_set$ref_snp[i] & finemap$feature == final_set$feature[i] &
					 #finemap$X.log_gwas_pval > 5 & finemap$X.log_eqtl_pval > 5,]
					 finemap$X.log_gwas_pval > 5 & finemap$X.log_eqtl_pval > 7.3,]

			if (dim(result)[1] == 0)
			{
				my_matrix[j,k] = NA
			}
			else
			{
				my_matrix[j,k] = result$clpp_mod[1]
			}
			if (dim(result_trunc)[1] == 0)
			{
				trunc_matrix[j,k] = NA
			}
			else
			{
				trunc_matrix[j,k] = result_trunc$clpp_mod[1]
			}
		}
	}

	# Figure out which rsid this is, from the original rsid file.

	# Make heatmaps for this matrix
	# First though, print it
	print(final_set[i,])
	print(my_matrix)

	plot_title = paste(final_set$rsid[i], final_set$ref_snp[i], final_set$feature[i], final_set$hgnc[i], sep="\n")

	# TODO:
	# Also consider changing color thresholds to avoid making them misleading...
	# to make them correspond more to our intuitions.

	# Plot heatmap
	heat = melt(my_matrix)
	heat$plot_vals = pmin(0.5, heat$value)
	heat$N = 0

	# Get the number of SNPs actually tested at this locus;
	# will be plotted within the heatmap tiles
	for (j in 1:dim(heat)[1])
	{
		item = finemap[finemap$eqtl_file == heat$X1[j] & finemap$gwas_trait == heat$X2[j] &
			       finemap$ref_snp == final_set$ref_snp[i] & finemap$feature == final_set$feature[i],]

		if (dim(item)[1] != 0)
		{
			heat$N[j] = item$n_snps
		}
	}	


	heat$text = paste(round(heat$value,3), "\n(n=", heat$N, ")", sep="")
	g = ggplot(data = heat, aes(x = X1, y = X2)) +
			geom_tile(aes(fill = plot_vals)) +
			geom_text(aes(label = text)) +
			scale_fill_gradient2(low = "white", high = 'orangered4', limits=c(0, 0.5)) +
			theme(axis.text.y=element_text(size=10),
			axis.text.x = element_text(angle = 90, hjust = 1)) +
	                labs(x = "eQTL Tissue") +
	                labs(y = "GWAS Trait") + 
			ggtitle(plot_title)
	
	ggsave(paste("/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/plots/clpp", final_set$ref_snp[i], "_", final_set$feature[i], "mod.png", sep=""), width = 8, height = 6, units = "in", dpi = 300, limitsize=FALSE)

	# Plot heatmap
	heat = melt(trunc_matrix)
	heat$plot_vals = pmin(0.5, heat$value)
	heat$N = 0

	# Get the number of SNPs actually tested at this locus
	for (j in 1:dim(heat)[1])
	{
		item = finemap[finemap$eqtl_file == heat$X1[j] & finemap$gwas_trait == heat$X2[j] &
			       finemap$ref_snp == final_set$ref_snp[i] & finemap$feature == final_set$feature[i],]

		if (dim(item)[1] != 0)
		{
			heat$N[j] = item$n_snps
		}
	}
	heat$text = paste(round(heat$value,3), "\n(n=", heat$N, ")", sep="")
	g = ggplot(data = heat, aes(x = X1, y = X2)) +
			geom_tile(aes(fill = plot_vals)) +
			geom_text(aes(label = text)) +
			scale_fill_gradient2(low = "white", high = 'orangered4', limits=c(0, 0.5)) +
			theme(axis.text.y=element_text(size=10),
			axis.text.x = element_text(angle = 90, hjust = 1)) +
	                labs(x = "eQTL Tissue") +
	                labs(y = "GWAS Trait") +
			ggtitle(plot_title)

	ggsave(paste("/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/plots/clpp", final_set$ref_snp[i], "_", final_set$feature[i], "_mod_cut-off.png", sep=""), width = 8, height = 6, units = "in", dpi = 300, limitsize=FALSE)

}

# NOTE: Should also collect LocusCompare-style plots on the side,
# to get a better idea of what's going on.

# Produce a list of genes for prioritization in follow-up studies, ranked by CLPP_mod score.
# Order by max CLPP score
# (Note that duplicates have already been removed earlier -- so this is how many hits we've got).
priority = final_set[rev(order(final_set$clpp_mod)),]

# Number of hits
print(dim(priority)[1])
# Number of SNPs that had at least one hit
print(length(unique(priority$rsid)))
# Contrast with total number of SNPs
print(length(unique(finemap$rsid)))

# Print a table showing SNP, gene, CLPP score in top context, coloc rank among all genes associated with that SNP.
# Other info can be found from the heatmaps

# We have most of this info already, but we need the column showing rank of this gene among all colocalized
# genes at this locus
priority$rank_at_locus = sapply(1:dim(priority)[1], function(x) 
	   	{
			return(sum(priority$rsid[1:x] == priority$rsid[x]))
		}
	   )

priority$hgnc = as.character(priority$hgnc)
priority[which(priority$hgnc == ""),]$hgnc = "missing"

# Then can plot the total number of colocalizations found for each gene and include that in the presentation too

last_only = priority[!rev(duplicated(rev(priority$rsid))),]	# Get last appearance of each rsid in the list
hist(last_only$rank_at_locus, breaks=0:max(last_only$rank_at_locus))

#  We'd also like to know, out of all genes tested in the region, how close is the gene compared to other nearby ones?
# We obtained this file from BioMart
gene_locs = read.csv("/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/ensembl_genes_with_coordinates.txt", header=TRUE, sep="\t")

priority$ensembl = as.character(priority$ensembl)
gene_locs$Gene.stable.ID = as.character(gene_locs$Gene.stable.ID)
gene_locs$Chromosome.scaffold.name = as.numeric(as.character(gene_locs$Chromosome.scaffold.name))
gene_locs = gene_locs[!is.na(gene_locs$Chromosome.scaffold.name),]

# Proximity score = proximity rank of the gene to our SNP, compared
# with all genes tested
priority$proximity = sapply(1:dim(priority)[1], function(x)
       {
		coloc_sub = priority[x,]
		coloc_chrom = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][1])
		coloc_pos = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][2])

       		# Get every nearby gene
       		gene_sub = gene_locs[gene_locs$Chromosome.scaffold.name == coloc_chrom,]
		# Subset down to only the set of genes that we actually tested for colocalization
	        gene_sub = gene_sub[gene_sub$Gene.stable.ID %in% finemap$ensembl,] 

		uniq = gene_sub %>% group_by(Gene.stable.ID) %>% summarize(start = min(c(Gene.start..bp., Gene.end..bp.)), end = max(c(Gene.start..bp., Gene.end..bp.)))

		# For all nearby genes, determine order of proximity to our SNP of interest

		dist = sapply(1:dim(uniq)[1], function(y)
		       {
				this_gene = uniq[y,]
		       		if (this_gene$start > coloc_pos)
				{
					return(this_gene$start - coloc_pos)
				}
				else if (this_gene$end < coloc_pos)
				{
					return(coloc_pos - this_gene$end)
				}
				else
				{
					return(0)
				}
		       }
		)

		uniq = uniq[order(dist),]

		return(which(uniq$Gene.stable.ID == coloc_sub$ensembl))
       })

# Get gene distance from TS
priority$tss_distance = sapply(1:dim(priority)[1], function(x)
       {
		coloc_sub = priority[x,]
		coloc_chrom = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][1])
		coloc_pos = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][2])

		my_gene = gene_locs[gene_locs$Gene.stable.ID == coloc_sub$ensembl,]
		return(min(abs(coloc_pos - my_gene$Transcription.start.site..TSS.)))
       })


# Could ask this question:
# If we consider all the genes that are colocalized with our SNP, what is the proximity
# rank of the BEST one? (This way, we're not penalizing genes that happen to just have a slightly
# lower colocalization score with the nearby gene than with the distant one).

min_proximity = priority %>% group_by(rsid) %>% summarize(min_proximity = min(proximity))
hist(min_proximity$min_proximity, breaks=0:max(min_proximity$min_proximity))

# Get total number of genes tested at this locus
genes_tested = finemap %>% group_by(rsid) %>% summarize(genes_tested=length(unique(as.character(ensembl))))
priority = merge(priority, genes_tested)

# Get total number of colocalizations found at this locus
last_only$num_colocs = last_only$rank_at_locus
priority = merge(priority, last_only[c("rsid", "num_colocs")], by=c("rsid"))

# Get ALL hits at this locus (regardless of whether others were found in same tissue
extra_info = priority[c("rsid", "ensembl", "rank_at_locus", "proximity", "tss_distance", "num_colocs", "genes_tested")]
priority = merge(pre_final_set, extra_info, by=c("rsid", "ensembl"))

num_gwas_colocalized = priority %>% group_by(rsid, ensembl) %>% summarize(num_gwas_colocalized = length(unique(base_gwas_file)))
priority = merge(priority, num_gwas_colocalized, by=c("rsid", "ensembl"))

num_eqtl_colocalized = priority %>% group_by(rsid, ensembl) %>% summarize(num_tissues_colocalized = length(unique(eqtl_file)))
priority = merge(priority, num_eqtl_colocalized, by=c("rsid", "ensembl"))

num_gwas_by_eqtl_colocalized = priority %>% group_by(rsid, ensembl) %>% summarize(num_gwas_by_eqtl_colocalized = length(unique(paste(base_gwas_file, eqtl_file))))
priority = merge(priority, num_gwas_by_eqtl_colocalized, by=c("rsid", "ensembl"))

priority = priority[rev(order(priority$clpp_mod)),]
write.table(priority, file="/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/results/prioritized_gene_rankings_no_clpp_mod_cutoff.tsv", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
