# Author: Mike Gloudemans
# Author: Jeremy Tien

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
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import csv
from sklearn.preprocessing import label_binarize
from scipy import interp
from sklearn.metrics import accuracy_score


def run_ensemble(locus):
    
    # TODO for Jeremy:
    # Compute ensemble score using scores passed in from other functions.
    # The scores from the other functions will be stored in locus.scores,
    # which is an object of class StoreContainer and should be the only information
    # you need for this I think. Scores that have not been
    # successfully computed (or just weren't requested) will have a value of NoneType.
    #
    
    ###################
    ## PREPROCESSING ##
    ###################

    # Load matrix of simulated data for training. Consider making this hard-coded in.
    matrix_save_path = locus.settings["methods"]["ensemble"]["training_data"] # current matrix contains COLOC, RTC, FINEMAP_STANDARD, FINEMAP_MODIFIED, CAVIARBF, BASELINE, SMARTBASELINE, SMR, SMR_HEIDIADJUSTED, GSMR, TWAS 
    matrix = pd.read_csv(matrix_save_path, sep="\t")

    # remove any submethods based on input
    methods_used = []
    if isnumeric(locus.scores.coloc_h4):
	methods_used.append("COLOC")
    if isnumeric(locus.scores.rtc_neg_log_pval): 
        methods_used.append("RTC")
    if isnumeric(locus.scores.finemap_clpp): 
        methods_used.append("FINEMAP_STANDARD")
    if isnumeric(locus.scores.finemap_clpp_mod): 
        methods_used.append("FINEMAP_MODIFIED")
    # CAVIARBF not yet implemented
    if isnumeric(locus.scores.baseline_neg_log_pval): 
        methods_used.append("BASELINE")
    if isnumeric(locus.scores.smart_baseline_neg_log_pval): 
        methods_used.append("SMARTBASELINE")
    if isnumeric(locus.scores.smr_neg_log_pval): 
        methods_used.append("SMR")
    # SMR_HEIDIADJUSTED not yet implemented
    if isnumeric(locus.scores.gsmr_neg_log_pval): 
        methods_used.append("GSMR")
    if isnumeric(locus.scores.twas_neg_log_pval): 
        methods_used.append("TWAS")
    
    matrix = matrix[methods_used]

    # use all data in matrix for training (no need for train test split)
    trainX = matrix.drop('colocalization_status', 1)
    trainY = matrix['colocalization_status']

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
