data = read.table("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/scz2_snp_results_txt_finemap_clpp_status.txt")

# Find the percentage of sites in each tissue whose coloc
# score exceeds a certain threshold.
clpp_threshold = 0.01
colocs = rep(0, length(levels(data[,2])))
png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/scz2_ecdf.png", units="in", width=5, height=5, res=360)
	plot(1, type="n", xlim=c(-2.75, -0.5), ylim=c(0, 1), ylab="eCDF of top 1% of tested sites in a tissue", xlab="CLPP score")
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
		print(10^top_one_percent)
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
	}
dev.off()

# Do a similar comparison, just with brain and non-brain
png("/users/mgloud/projects/brain_gwas/output/scz2_snp_results_txt/1.0e-07/scz2_brain_nonbrain_ecdf.png", units="in", width=5, height=5, res=360)

plot(1, type="n", xlim=c(-2.75, -0.5), ylim=c(0, 1), ylab="eCDF of top 1% of tested sites in a tissue", xlab="CLPP score")

brain_colocs = data[grep("Brain", data[,2]),]
brain_clpp = brain_colocs[,6]
percent_brain_colocs = sum(brain_clpp > clpp_threshold) / length(brain_clpp) * 100
brain_top_one_percent = log10(brain_clpp[order(brain_clpp)][(round(length(brain_clpp)*99/100)):length(brain_clpp)])
lines(brain_top_one_percent,(1:length(brain_top_one_percent))/length(brain_top_one_percent),type="s")

nonbrain_colocs = data[-grep("Brain", data[,2]),]
nonbrain_clpp = nonbrain_colocs[,6]
percent_nonbrain_colocs = sum(nonbrain_clpp > clpp_threshold) / length(nonbrain_clpp) * 100
nonbrain_top_one_percent = log10(nonbrain_clpp[order(nonbrain_clpp)][(round(length(nonbrain_clpp)*99/100)):length(nonbrain_clpp)])
lines(nonbrain_top_one_percent,(1:length(nonbrain_top_one_percent))/length(nonbrain_top_one_percent),col="red")

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
