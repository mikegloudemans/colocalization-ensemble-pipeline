#!/usr/bin/Rscript
# Date created: 11/27/2017
# Author: Mike Gloudemans

library(readr)
library(dplyr)

# Prep data
genotype = read_delim(gzfile("/srv/persistent/bliu2/rpe/data/genotype/asvcf/glucose_nodup/rpe.imputed.chr12.all_filters.vcf.new.gz"), comment="##", delim="\t")
genotype1 = genotype[grep("_56115585_", genotype$ID),]
genotype2 = genotype[grep("_56115778_", genotype$ID),]

genotype1 = t(genotype1)
g1_table = data.frame(id=rownames(genotype1)[10:dim(genotype1)[1]], variant=substr(genotype1[10:dim(genotype1)[1],1],1,3))
genotype2 = t(genotype2)
g2_table = data.frame(id=rownames(genotype2)[10:dim(genotype2)[1]], variant=substr(genotype2[10:dim(genotype2)[1],1],1,3))

splicing = read_delim(gzfile("/srv/persistent/bliu2/rpe/data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr12.gz"), delim="\t")
splice = splicing[grep("12:56115278:56117670:clu_1202", splicing$ID),]
splice = t(splice)
splice_table = data.frame(id=rownames(splice)[5:dim(splice)[1]], splicing=as.numeric(splice[5:dim(splice)[1],1]))

expression = read_delim("/srv/persistent/bliu2/rpe/data/rnaseq/count/merged/rpe.gene_count", delim="\t")
rdh5 = expression[grep("ENSG00000135437", expression$gene_id),]
glucose_rdh5 = rdh5[colnames(rdh5)[grep("glucose", colnames(rdh5))]]
galactose_rdh5 = rdh5[colnames(rdh5)[grep("galactose", colnames(rdh5))]]
colnames(glucose_rdh5) = gsub(".glucose", "", colnames(glucose_rdh5))
colnames(galactose_rdh5) = gsub(".galactose", "", colnames(galactose_rdh5))

glucose_rdh5 = t(glucose_rdh5)
glucose_rdh5_table = data.frame(id=rownames(glucose_rdh5), expression=glucose_rdh5[,1])
galactose_rdh5 = t(galactose_rdh5)
galactose_rdh5_table = data.frame(id=rownames(galactose_rdh5), expression=galactose_rdh5[,1])

# The two genotype tables are the same, so just ditch the second one
print(inner_join(g1_table, g2_table, by="id"))

all_data = inner_join(inner_join(g1_table, splice_table), glucose_rdh5_table)
all_data$variant = all_data$variant != "0|0"
all_data$expression  = log(all_data$expression)

# Analysis
png("plots/RDH5_sqtl_12_56115585.png", units="in", height=5, width=5, res=360)
	boxplot(splicing ~ variant, data=all_data)
dev.off()

boxplot(expression ~ variant, data=all_data)
plot(all_data$expression, all_data$splicing)

t.test(all_data[all_data$variant,]$expression, all_data[!all_data$variant,]$expression)

