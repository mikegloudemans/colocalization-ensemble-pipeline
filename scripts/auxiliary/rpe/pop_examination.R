d = read_table("1kg_alleles.txt")

ids = read_table("ATACSeqSampleIDs.txt", col_names=FALSE)
samps = sapply(ids$X1, function(x) {substring(x,1,7)})
names(samps) = 1:length(samps)

# None of the African subset (samples we have for African
# genomes project) is anything but homozygous reference
african_subset = d[,which(colnames(d) %in% samps)]
as.matrix(african_subset)

homo_wt = d[,which(d[1,] == "1|1" | d[2,] == "1|1")]
homo_ref = d[,which(d[1,] == "0|0" | d[2,] == "0|0")]

pops = read_csv("20130606_sample_info.csv")
homo_wt_pops = pops[match(names(homo_wt), pops$Sample),]

table(homo_wt_pops$Population)
table(pops$Population)

pop_membs = pops[match(colnames(d), pops$Sample),]$Population
colnames(d)[10:length(colnames(d))] = paste(colnames(d)[10:length(colnames(d))], pop_membs[10:length(pop_membs)], sep="_")

write_delim(homo_wt_pops, path="homo_wt_individuals.csv", delim=",")


