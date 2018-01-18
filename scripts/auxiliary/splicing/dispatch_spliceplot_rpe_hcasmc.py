import run_splice_plot as sp

hcasmc_results_file = "/users/mgloud/projects/brain_gwas/output/2017-12-13_11-47-47_hcasmc_splicing/UKBB_GWAS1KG_EXOME_CAD_SOFT_META_PublicRelease_300517_txt_gz_finemap_clpp_status.txt"
hcasmc_splice_file = "/users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/sqtls/hcasmc_sqtls.txt.gz"
hcasmc_map_file = "/users/mgloud/projects/brain_gwas/data/spliceplot/hcasmc_map_file.txt"
hcasmc_vcf_file = "/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/asvcf/phased_and_imputed.chr{0}.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz"

sp.generate_splice_plots(hcasmc_results_file, hcasmc_splice_file, hcasmc_vcf_file, hcasmc_map_file, threshold = 0.003, junction_delim="_")
#sp.generate_splice_plots(hcasmc_results_file, hcasmc_splice_file, hcasmc_vcf_file, hcasmc_map_file, threshold = 0.005, junction_delim="_")

rpe_results_file = "/users/mgloud/projects/brain_gwas/output/2017-11-14_07-52-48_rpe_splicing/Fritsche_sorted_txt_gz_finemap_clpp_status.txt"
rpe_splice_file = "/users/mgloud/projects/brain_gwas/data/eqtls/rpe/sqtls/rpe_sqtls.txt.gz"
rpe_map_file = "/users/mgloud/projects/brain_gwas/data/spliceplot/rpe_glucose_map_file.txt"
rpe_vcf_file = "/srv/persistent/bliu2/rpe/data/genotype/asvcf/glucose_nodup/rpe.imputed.chr{0}.all_filters.vcf.new.gz"

sp.generate_splice_plots(rpe_results_file, rpe_splice_file, rpe_vcf_file, rpe_map_file, threshold = 0.005)


