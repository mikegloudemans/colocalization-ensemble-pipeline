# this script is not properly commented, it just makes SQLite DBs from the intermediate files required for this step
# for real data, this entire section can be skipped if previously created DBs are used, which we would use on real data
# alternately, this script can still work for real data, but genotypes (and other required files such as expression levels and covariates) would need to be provided

library(dplyr)
library(RSQLite)
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
tissue <- argv[1] # needs to be sim{sim_num}
chrom <- argv[2]
gene_annot_file <- argv[3] # /users/raoa/coloc_comparison/tmp.annot286.gtf
select_signif <- as.logical(argv[4]) # whether final db is filtered by zscore_pval and rho_avg, in many cases in the simulation, there will be no genes remaining, so this is set to False in the main script
working_directory <- argv[5]

# dir = working_directory %&% 'PredictDB_Pipeline_GTEx_v7/model_training/scripts'
# setwd(dir)

driver <- dbDriver('SQLite')
gene_annot <- read.table(gene_annot_file, header = T, stringsAsFactors = F)

model_summaries <- read.table('../summary/' %&% tissue %&% '_nested_cv_chr' %&% as.character(chrom) %&% '_model_summaries.txt', header = T, stringsAsFactors = F)
tiss_summary <- read.table('../summary/' %&% tissue %&% '_nested_cv_chr' %&% as.character(chrom) %&% '_tiss_chr_summary.txt', header = T, stringsAsFactors = F)
  
n_samples <- tiss_summary$n_samples

model_summaries <- rename(model_summaries, gene = gene_id)
if (file.exists('../dbs/' %&% tissue %&% '_tw0.5.db')) {system('rm ../dbs/' %&% tissue %&% '_tw0.5.db')}
conn <- dbConnect(drv = driver, '../dbs/' %&% tissue %&% '_tw0.5.db')
dbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)
dbGetQuery(conn, "CREATE INDEX gene_model_summary ON model_summaries (gene)")

weights <- read.table('../weights/' %&% tissue %&% '_nested_cv_chr' %&% as.character(chrom) %&%'_weights.txt', header = T, stringsAsFactors = F)
weights <- rename(weights, gene = gene_id)

dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbGetQuery(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbGetQuery(conn, "CREATE INDEX weights_gene ON weights (gene)")
dbGetQuery(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

sample_info <- data.frame(n_samples = n_samples, population = 'simulated_european', tissue = tissue)
dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)

construction <- tiss_summary %>% select(chrom, cv_seed) %>% rename(chromosome = chrom)
dbWriteTable(conn, 'construction', construction, overwrite = TRUE)

# select significant? and create final db

if (file.exists('../dbs/' %&% tissue %&% '_tw0.5_filtered.db')) {system('rm ../dbs/' %&% tissue %&% '_tw0.5_filtered.db')}
out_conn <- dbConnect(driver, '../dbs/' %&% tissue %&% '_tw0.5_filtered.db')
if (select_signif) {
	model_summaries <- dbGetQuery(conn, 'select * from model_summaries where zscore_pval < 0.05 and rho_avg > 0.1') 
} else {
	model_summaries <- dbGetQuery(conn, 'select * from model_summaries')
}
model_summaries <- model_summaries %>% rename(pred.perf.R2 = rho_avg_squared, genename = gene_name, pred.perf.pval = zscore_pval, n.snps.in.model = n_snps_in_model)
model_summaries$pred.perf.qval <- NA

dbWriteTable(out_conn, 'extra', model_summaries)
construction <- dbGetQuery(conn, 'select * from construction')
dbWriteTable(out_conn, 'construction', construction)
sample_info <- dbGetQuery(conn, 'select * from sample_info')
dbWriteTable(out_conn, 'sample_info', sample_info)
weights <- dbGetQuery(conn, 'select * from weights')
weights <- weights %>% filter(gene %in% model_summaries$gene) %>% rename(eff_allele = alt, ref_allele = ref, weight = beta)
dbWriteTable(out_conn, 'weights', weights)
dbGetQuery(out_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbGetQuery(out_conn, "CREATE INDEX weights_gene ON weights (gene)")
dbGetQuery(out_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
dbGetQuery(out_conn, "CREATE INDEX gene_model_summary ON extra (gene)")

