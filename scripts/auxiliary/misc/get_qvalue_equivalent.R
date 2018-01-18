# Author: Mike Gloudemans
# Date created: 12/7/2017
# Get the p-value corresponding to a q-value of 0.05 in eCAVIAR input

require(readr)
require(qvalue)

d = read_delim(gzfile("UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz"), delim="\t")
d$qvalue = qvalue(d$pvalue)$qvalues
sub = d[d$qvalue < 0.05,]
print(sub[rev(order(sub$qvalue)),][c("pvalue", "qvalue")])
