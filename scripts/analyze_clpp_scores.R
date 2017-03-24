data = read.table("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/fixed_threshold/1.0e-06_by_1.0e-05/scz2_snp_results_txt_clpp_status.txt")
#data = read.table("/users/mgloud/projects/brain_gwas/output/pgc_mdd_full_2012-04_hg19_txt/fixed_threshold/1.0e-06_by_5.0e-04/pgc_mdd_full_2012-04_hg19_txt_clpp_status.txt")
#data = read.table("/users/mgloud/projects/brain_gwas/output/anxiety_meta_full_fs_tbl/fixed_threshold/1.0e-06_by_5.0e-05/anxiety_meta_full_fs_tbl_clpp_status.txt")

# Find the percentage of sites in each tissue whose coloc
# score exceeds a certain threshold.
clpp_threshold = 0.01
colocs = rep(0, length(levels(data[,2])))
png("/users/mgloud/projects/brain_gwas/output/scz2_ecdf.png", units="in", width=5, height=5, res=360)
	plot(1, type="n", xlim=c(-2.25, -0.5), ylim=c(0, 1), xlab="eCDF of top 1% of tested sites in a tissue", ylab="CLPP score")
	for (i in 1:length(levels(data[,2])))
	{
		tissue = levels(data[,2])[i]
		subset = data[data[,2] == tissue,]
		clpp = subset[,6]
		percent_colocs = sum(clpp > clpp_threshold) / length(clpp) * 100
		colocs[i] = percent_colocs

		# Plot ecdf for top 1% of coloc scores
		top_one_percent = log10(clpp[order(clpp)][(round(length(clpp)*99/100)):length(clpp)])
		lines(top_one_percent,(1:length(top_one_percent))/length(top_one_percent),type="s")
		print(top_one_percent)
		print(exp(top_one_percent))
		#readline(prompt=tissue)
	}
	coloc_table = cbind(levels(data[,2]),colocs)
	coloc_table = coloc_table[rev(order(coloc_table[,2])),]
	print(coloc_table)

	# Just for simplicity of visualization, plot brain tissues on top in a different color
	colors = rainbow(length(grep("Brain", levels(data[,2]))))
	brain_index = 0
	for (i in 1:length(levels(data[,2])))
	{
		tissue = levels(data[,2])[i]
		if (length(grep("Brain", tissue)) != 1)
		{
			next;
		} 
		brain_index = brain_index + 1
		subset = data[data[,2] == tissue,]
		clpp = subset[,6]
		
		# Plot ecdf for top 1% of coloc scores
		top_one_percent = log10(clpp[order(clpp)][round((length(clpp)*99/100)):length(clpp)])
		lines(top_one_percent,(1:length(top_one_percent))/length(top_one_percent),type="s", col="black", lwd=6)
		lines(top_one_percent,(1:length(top_one_percent))/length(top_one_percent),type="s", col=colors[brain_index], lwd=3)
		#readline(prompt=tissue)
	}
dev.off()

# Observations: Not much evidence that brain tissues show greater % colocalization than
# non-brain tissues, with ~130 GWAS SNPs analyzed. Some of top tissues are indeed brain tissues,
# but so are some of the bottom ones.
#
# For lower thresholds (such as 0.001 required CLPP score), brain tissues do appear towards
# the top of the list, BUT so do other tissues of low sample size (Uterus, Vagina, Ovary).
# There seems to be a sample size dependency going on once we get down to lower levels.
#
# My inclination is that we should only trust colocalization scores if there's quite high confidence
# in the score. Even then, we need to take it with a grain of salt because it still may be somewhat
# susceptible to sample size effects.
# 
# One other potential issue with sample size effects: Suppose SNPs are profiled more densely for eQTL
# in a non-brain tissue. Then, all other things equal, if there is no signal, the brain tissue may
# receive a slightly higher coloc score simply due to the fact that there are more tied SNPs in the
# non-brain tissue. (I'll have to think about this a bit more though - this may not be entirely accurate).
#
# I should collect my thoughts at this and present soon at a lab meeting.
#
