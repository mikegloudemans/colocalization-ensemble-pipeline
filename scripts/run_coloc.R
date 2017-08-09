# TODO: Along with the Python dispatcher,
# come up with a standard input format for input files.



# Get input and output locations as command-line arguments

infile = sys.argv[1]
outfile = sys.argv[2]

# Run coloc
data = read.table(infile, header=TRUE)

# TODO: Fix parameter settings in here. This is just for an initial run-through.
# GWAS data
dataset1 = list(pvalues=data$pvalue_gwas, MAF=data$allele_freq, s=(data$Ncases[1] / (data$Ncases[1] + data$Ncontrols[1])), snp=data$rsid, N=data$Ncases+data$Ncontrols, type="cc")
# eQTL data
dataset2 = list(pvalues=data$pvalue_eqtl, MAF=data$allele_freq, snp=data$rsid, N=50, type="quant")

results = coloc.abf(dataset1, dataset2)

# TODO: Output files to desired location.
# Will be parsed further in Python.

# (Actually though...will probably be easier to just output to the command line, for reading into Python. All I have to do is read H4.)

write.table(results, file=outfile)
