# Mike Gloudemans
# Date Created: 4/18/2018

data = read.table("/users/mgloud/projects/brain_gwas/output/2018-04-18_09-42-40_insulin_for_eriks_grant/WHRadjBMI_GIANT_Mixed_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt", header = TRUE)
data = data[data$X.log_gwas_pval > 5 & data$X.log_eqtl_pval > 5,]

data = data[data$eqtl_file == "Muscle_Skeletal_allpairs_txt_gz",]


plot(data$clpp, data$clpp_mod)
plot(log10(data$clpp), log10(data$clpp_mod))

data_coloc = data[data$clpp > 0.01 | data$clpp_mod > 0.1,]
plot(data_coloc$clpp, data_coloc$clpp_mod)
plot(log10(data_coloc$clpp), log10(data_coloc$clpp_mod))

