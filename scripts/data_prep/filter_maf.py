# Author: Mike Gloudemans
# Filter variants from VCF if MAF < 0.005
# Yes, VCFtools does this too, but just slower

import sys
import gzip

input = sys.argv[1]
output = sys.argv[2]

with open(output, "w") as w:
    with gzip.open(input, "rb") as f:
        for line in f:
            if line.startswith("#"):
                w.write(line)
            else:
                af = line.split(";AF=")[1].split(";")[0]
                if "," in af:
                    continue
                af = float(af)
                if af < 0.005 or 1-af < 0.005:
                    continue
                else:
                    w.write(line)
