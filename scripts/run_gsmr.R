#!/usr/bin/env Rscript
# Author: Abhiram Rao
# Author: Mike Gloudemans
# requires GCTA http://cnsgenomics.com/software/gcta/gcta_1.92.1beta6.zip

suppressMessages(require(tidyverse))
suppressMessages(require(gsmr))
suppressMessages(require(survey))
suppressMessages(require(data.table))

sink(file="/dev/null")

# We will assume that the wrapper script (in Python)
# has written out the combined data frame to the specified files
# and has specified the other relevant parameters

# sim number as argument
args <- commandArgs(trailingOnly = TRUE)


tmp_base = args[1]
combined_data = read_csv(paste0(tmp_base, ".csv"))

# options required - config?
eqtl_data = '/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/eqtl'
gwas_data = '/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/gwas'
ld_reference = '/users/raoa/coloc_comparison/ld_reference'
gsmr_result_location = '/users/raoa/coloc_comparison/gsmr_results'

# Parameters for running: may add options to modify these later.
gcta_location = '/users/raoa/coloc_comparison/gcta_1.92.1beta6/gcta64'
standardize = FALSE #If the risk factor was not standardised in GWAS, the effect sizes can be scaled, that this process requires allele frequencies, z-statistics and sample size. After scaling, bzx is interpreted as the per-allele effect of a SNP on the exposure in standard deviation units
n_ref = 2504    # Sample size of the reference sample
#gwas_thresh = 5e-5    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
gwas_thresh = 5e-3    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
single_snp_heidi_thresh = 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
multi_snp_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
nsnps_thresh = 2   # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag = FALSE    # flag for HEIDI-outlier analysis
#ld_r2_thresh = 0.2    # LD r2 threshold to remove SNPs in high LD
ld_r2_thresh = 0.05    # LD r2 threshold to remove SNPs in high LD
ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
gsmr_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development 

chrom = unique(combined_data$chr_eqtl)[1]

# NOTE: This part could be made faster by never splitting the data to eqtl and gwas in the first place
# subset to eqtl and gwas summary data
eqtl = combined_data[c("rsid", "ref_eqtl", "alt_eqtl", "effect_af_eqtl", "beta_eqtl", "se_eqtl", "pvalue_eqtl", "N")]
#gene = unique(eqtl$feature)[1]
gwas = combined_data[c("rsid", "beta_gwas", "se_gwas", "pvalue_gwas", "n_cases", "n_controls")]
gwas = mutate(gwas, N = n_cases + n_controls)
gwas = combined_data[c("rsid", "beta_gwas", "se_gwas", "pvalue_gwas", "N")]
colnames(eqtl) = c("SNP", "a2", "a1", "a1_freq", "bzx", "bzx_se", "bzx_pval", "bzx_n")
colnames(gwas) = c("SNP", "bzy", "bzy_se", "bzy_pval", "bzy_n")

gsmr_df = merge(eqtl, gwas, by="SNP")
if (dim(gsmr_df)[1] == 0) {stop("No variants shared between eQTL/GWAS")}

# calculate LD correlation matrix using GCTA
write.table(gsmr_df[,c(1,2)], sprintf(paste0(tmp_base, ".allele")), col.names=F, row.names=F, quote=F)
system(sprintf("%s --bfile %s/EUR1KG_chr%s_hg19_rsid --extract %s.allele --update-ref-allele %s.allele --recode --out %s.ld_mat > /dev/null", gcta_location, ld_reference, chrom, tmp_base, tmp_base, tmp_base))

snp_coeff_id = scan(sprintf("%s.ld_mat.xmat.gz", tmp_base), what="", nlines=1)
snp_coeff = read.table(sprintf("%s.ld_mat.xmat.gz", tmp_base), header=F, skip=2)

# Match the SNP genotype data with the summary data
snp_id = Reduce(intersect, list(gsmr_df$SNP, snp_coeff_id))
gsmr_df = gsmr_df[match(snp_id, gsmr_df$SNP),]
snp_order = match(snp_id, snp_coeff_id)
snp_coeff_id = snp_coeff_id[snp_order]
snp_coeff = snp_coeff[, snp_order]

# Calculate the LD correlation matrix
ldrho = cor(snp_coeff)

# Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
colnames(ldrho) = rownames(ldrho) = snp_coeff_id

if (standardize) {
	snpfreq = gsmr_df$a1_freq             # allele frequencies of the SNPs
	bzx = gsmr_df$bzx     # effects of the instruments on risk factor
	bzx_se = gsmr_df$bzx_se       # standard errors of bzx
	bzx_n = gsmr_df$bzx_n          # GWAS sample size for the risk factor
	std_zx = std_effect(snpfreq, bzx, bzx_se, bzx_n)    # perform standardisation
	gsmr_df$std_bzx = std_zx$b    # standardized bzx
	gsmr_df$std_bzx_se = std_zx$se    # standardized bzx_se
}

## run GSMR
# saves a file with just p-value if GSMR runs successfully, there is no file generated otherwise

pval = 1

tryCatch(
	 {
		gsmr_results = suppressMessages(gsmr(gsmr_df$bzx, gsmr_df$bzx_se, gsmr_df$bzx_pval, gsmr_df$bzy, gsmr_df$bzy_se, gsmr_df$bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr_beta))
		pval = gsmr_results$bxy_pval
	 },
	 error = function(err){})

sink()

cat(pval)

## NOTE
# bzx is SNP effect on risk factor, in this case gene expression
# bzy is the SNP effect on the disease, from GWAS
