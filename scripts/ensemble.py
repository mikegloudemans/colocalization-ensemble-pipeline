# Author: Mike Gloudemans

import subprocess
from scipy import stats
from shutil import copyfile
import sys
if sys.version_info[0] < 3: 
   from StringIO import StringIO
else:
   from io import StringIO
import pandas as pd
import math

def run_baseline(locus):
    
    # TODO for Jeremy:
    # Compute ensemble score using scores passed in from other functions.
    # The scores from the other functions will be stored in locus.scores,
    # which is an object of class StoreContainer and should be the only information
    # you need for this I think. Scores that have not been
    # successfully computed (or just weren't requested) will have a value of NoneType.
    #

    # TODO: Set this value equal to the ensemble score once you've computed it, 
    # and then it will be written to the output file.
    ensemble_score = 0

    # TODO for Jeremy: Write the names of all the methods / scores that were included (not NoneType) in this file so we can
    # reference it later
    with open("{0}/ensemble_methods_used.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as w:
        pass

    with open("{0}/{1}_ensemble_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, locus.data.shape[0], locus.gwas_suffix, ensemble_score))

    return ensemble_score
