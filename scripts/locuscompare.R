require(locuscomparer)
require(ggplot2)

args = commandArgs(trailingOnly=TRUE)

file1 = args[1]
file2 = args[2]
out_file = args[3]
title1 = args[4]
title2 = args[5]
vcf = args[6]
snp = args[7]

file1 = read.table(file1, header=TRUE)
file2 = read.table(file2, header=TRUE)

m = main(in_fn1 = file1, in_fn2 = file2, title=title1, title2=title2, snp=snp) # vcf_fn=vcf, snp=snp)
m

# Plot with low quality to reduce file size; re-plot later if needing a high-quality figure
ggsave(out_file, width=12, height=6, dpi=100)

