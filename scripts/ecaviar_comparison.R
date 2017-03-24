require(reshape2)

# Read in data
d = read.table("/users/mgloud/projects/brain_gwas/output/prefiltering_test/scz2_snp_results_txt/fixed_threshold/1.0e-05_by_1.0e-07/scz2_snp_results_txt_clpp_status_prefiltering_test.txt")
d$id = paste(d[,1], d[,2], d[,3], sep="_")

d_merged <- dcast(d, id~V7, value.var="V6")
png("/users/mgloud/projects/brain_gwas/output/ecaviar_comparison.png", units="in", width=5, height=5, res=360)
plot(d_merged[,2], d_merged[,3], xlab="eCAVIAR CLPP score, before filtering", ylab="eCAVIAR CLPP score, after filtering", xlim = c(0,0.02), ylim=c(0,0.02))
abline(a=0,b=1,col="red")
dev.off()

png("/users/mgloud/projects/brain_gwas/output/ecaviar_log_comparison.png", units="in", width=5, height=5, res=360)
plot(log(d_merged[,2]), log(d_merged[,3]), xlab="eCAVIAR log-CLPP score, before filtering", ylab="eCAVIAR log-CLPP score, after filtering")
abline(a=0,b=1,col="red")
dev.off()
