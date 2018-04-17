# Author: Mike Gloudemans
# Date created: 4/10/2018

require(reshape)
require(ggplot2)

#
# Load files
#
#files = dir("/users/mgloud/projects/brain_gwas/output/2018-04-10_11-24-04_insulin_for_eriks_grant")
files = dir("/users/mgloud/projects/brain_gwas/output/2018-04-13_14-26-25_insulin_for_eriks_grant")
files = files[grepl("finemap_clpp_status.txt", files)]

# Finemap results
finemap_results = list()
for (i in 1:length(files))
{
	#finemap_results[[i]] = read.table(gzfile(paste("/users/mgloud/projects/brain_gwas/output/2018-04-10_11-24-04_insulin_for_eriks_grant/", files[i], sep="")), header=TRUE)
	finemap_results[[i]] = read.table(gzfile(paste("/users/mgloud/projects/brain_gwas/output/2018-04-13_14-26-25_insulin_for_eriks_grant/", files[i], sep="")), header=TRUE)
}
finemap = do.call(rbind, finemap_results)
finemap$ensembl = sapply(as.character(finemap$feature), function(x) {strsplit(x, "\\.")[[1]][1]})
finemap$gwas_file = sapply(as.character(finemap$gwas_file), function(x) 
       {
	       s = strsplit(x, "_")[[1]]
	       return(paste(s[1:(length(s)-3)], collapse="_"))
       })

finemap$eqtl_file = sapply(as.character(finemap$eqtl_file), function(x)
       {
               s = strsplit(x, "_")[[1]]
               return(paste(s[1:(length(s)-3)], collapse="_"))
       })

# Load rsIDs
rsids = read.table("snps_with_hg19_positions.txt", header=FALSE)
rsids[,2] = gsub("chr", "", rsids[,2])
rsids$ref_snp = paste(rsids[,2], rsids[,3], sep="_")
rsids = rsids[,-c(2,3)]
colnames(rsids)[1] = "rsid"

# Get HGNC gene names
genes = read.table("mart_export.txt", header=TRUE, sep="\t")
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
dim(pass_pval_cutoffs)	# Total number of tests (gwas SNP - eqtl gene pairs) passing eQTL and GWAS pval thresholds
final_set = pass_pval_cutoffs[pass_pval_cutoffs$clpp > 0.1,]
dim(final_set)	# Total number of tests passing our colocalization criteria

# Remove duplicates from this list
final_set = final_set[!duplicated(paste(final_set$ref_snp, final_set$feature, sep="_")),]

#
# Then, for each SNP-gene pair passing the test, show full summary statistics across
# all eQTL files and GWAS
#
eqtl_files = unique(finemap$eqtl_file)
gwas_files = unique(finemap$gwas_file)

# For each SNP-gene pair that had at least one colocalization
for (i in 1:dim(final_set)[1])
{
	# Make a matrix to show colocalization across all tissues
	my_matrix = array(0, c(length(eqtl_files), length(gwas_files)))
	rownames(my_matrix) = eqtl_files
	colnames(my_matrix) = gwas_files
	trunc_matrix = my_matrix	# This matrix will only include the SNPs if they pass the p-value cutoffs

	# For each eQTL file
	for (j in 1:length(eqtl_files))
	{
		# For each GWAS file
		for (k in 1:length(gwas_files))
		{
			# Fill matrix accordingly, or fill with NA
			result = finemap[finemap$eqtl_file == eqtl_files[j] & finemap$gwas_file == gwas_files[k] & 
					 finemap$ref_snp == final_set$ref_snp[i] & finemap$feature == final_set$feature[i],]
			result_trunc = finemap[finemap$eqtl_file == eqtl_files[j] & finemap$gwas_file == gwas_files[k] & 
					 finemap$ref_snp == final_set$ref_snp[i] & finemap$feature == final_set$feature[i] &
					 finemap$X.log_gwas_pval > 5 & finemap$X.log_eqtl_pval > 5,]

			if (dim(result)[1] == 0)
			{
				my_matrix[j,k] = NA
			}
			else
			{
				my_matrix[j,k] = result$clpp[1]
			}
			if (dim(result_trunc)[1] == 0)
			{
				trunc_matrix[j,k] = NA
			}
			else
			{
				trunc_matrix[j,k] = result_trunc$clpp[1]
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
		item = finemap[finemap$eqtl_file == heat$X1[j] & finemap$gwas_file == heat$X2[j] &
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
	
	ggsave(paste("plots/clpp_mod", final_set$ref_snp[i], "_", final_set$feature[i], ".png", sep=""), width = 8, height = 6, units = "in", dpi = 300, limitsize=FALSE)

	# Plot heatmap
	heat = melt(trunc_matrix)
	heat$plot_vals = pmin(0.5, heat$value)
	heat$N = 0

	# Get the number of SNPs actually tested at this locus
	for (j in 1:dim(heat)[1])
	{
		item = finemap[finemap$eqtl_file == heat$X1[j] & finemap$gwas_file == heat$X2[j] &
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

	ggsave(paste("plots/clpp_mod", final_set$ref_snp[i], "_", final_set$feature[i], "_cut-off.png", sep=""), width = 8, height = 6, units = "in", dpi = 300, limitsize=FALSE)

}

# NOTE: Should also collect LocusCompare-style plots on the side,
# to get a better idea of what's going on.
