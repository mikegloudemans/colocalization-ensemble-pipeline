# Author: Mike Gloudemans

class SNP:
    def __init__(self, snp):
        if "chr" in str(snp[0]):
            self.chrom = int(snp[0][3:])
        else:
            self.chrom = int(snp[0])
        self.pos = int(snp[1])
        # we don't always need p-value
        if len(snp) > 2:
            self.pval = float(snp[2])
        else:
            self.pval = 0
