# Author: Mike Gloudemans

class SNP:
    def __init__(self, snp):
        if "chr" in str(snp[0]):
            self.chrom = int(snp[0][3:])
        else:
            self.chrom = int(snp[0])
        self.pos = int(snp[1])
        self.pval = float(snp[2])


