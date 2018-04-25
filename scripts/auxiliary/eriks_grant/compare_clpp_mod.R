# Mike Gloudemans
# Date Created: 4/18/2018

#
# Compare with just WHR in skeletal muscle
#

data = read.table("/users/mgloud/projects/brain_gwas/output/2018-04-18_09-42-40_insulin_for_eriks_grant/WHRadjBMI_GIANT_Mixed_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt", header = TRUE)
data = data[data$X.log_gwas_pval > 5 & data$X.log_eqtl_pval > 5,]

data = data[data$eqtl_file == "Muscle_Skeletal_allpairs_txt_gz",]


plot(data$clpp, data$clpp_mod)
plot(log10(data$clpp), log10(data$clpp_mod))

data_coloc = data[data$clpp > 0.01 | data$clpp_mod > 0.1,]
data = data[data$X.log_gwas_pval > 5 & data$X.log_eqtl_pval > 5,]
plot(data_coloc$clpp, data_coloc$clpp_mod)
plot(log10(data_coloc$clpp), log10(data_coloc$clpp_mod))

#
# Compare with multiple traits in multiple tissues
#
files = c("/users/mgloud/projects/brain_gwas/output/2018-04-19_07-49-06_insulin_for_eriks_grant/FastGlu_MAGIC_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/2018-04-19_07-49-06_insulin_for_eriks_grant/FastInsu_adjBMI_MAGIC_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/2018-04-19_07-49-06_insulin_for_eriks_grant/HDL_GCLC_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/2018-04-19_07-49-06_insulin_for_eriks_grant/TG_GCLC_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
		"/users/mgloud/projects/brain_gwas/output/2018-04-19_07-49-06_insulin_for_eriks_grant/WHRadjBMI_GIANT_Mixed_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt"
	  )

results = list()
for (file in files)
{
	results[[file]] = read.table(file, header=TRUE)
}

results = do.call(rbind, results)
rownames(results) = 1:dim(results)[1]

results = results[results$X.log_gwas_pval > 5 & results$X.log_eqtl_pval > 5,]
results_coloc = results[results$clpp > 0.01 | results$clpp_mod > 0.1,]


# Tests: how is CLPP score affected by number of SNPs tested?
# How is CLPP_mod score affected by number of SNPs tested?

plot(results$n_snps, results$clpp)

# In each bin of width 100 for n_snps, what percentage of tests have colocalization (score > 0.02)?
coloc_pcts = sapply(seq(0,1900,by=100), function(x)
       {
		subset = results[results$n_snps >= x & results$n_snps < x+100,]
       		hits = subset[subset$clpp > 0.01,]

		return(dim(hits)[1] / dim(subset)[1])
       })

lines(seq(0,1900,by=100) + 50, coloc_pcts, col="red", lwd=5)


plot(log2(results$n_snps), results$clpp)
lines(log2(seq(0,1900,by=100) + 50), coloc_pcts, col="red", lwd=5)

# Now do it with CLPP_mod
plot(results$n_snps, results$clpp_mod)

# In each bin of width 100 for n_snps, what percentage of tests have colocalization (score > 0.02)?
coloc_pcts = sapply(seq(0,1900,by=100), function(x)
       {
		subset = results[results$n_snps >= x & results$n_snps < x+100,]
       		hits = subset[subset$clpp_mod > 0.2,]

		return(dim(hits)[1] / dim(subset)[1])
       })

lines(seq(0,1900,by=100) + 50, coloc_pcts, col="red", lwd=5)


plot(log2(results$n_snps), results$clpp_mod)
lines(log2(seq(0,1900,by=100) + 50), coloc_pcts, col="red", lwd=5)

# Conclusion...I'm not sure if this is the perfect way to analyze it,
# but based on this preliminary analysis it doesn't seem obvious that the CLPP_mod
# score is less biased by number of SNPs included.

# Might also just be too small of a set of events to test this on though

