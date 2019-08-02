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
import numbers
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import csv
from sklearn.preprocessing import label_binarize
from scipy import interp
from sklearn.metrics import accuracy_score


def run_ensemble(locus):
    # TODO: rectify package versions and corresponding dependencies later
    # when updating software to latest.
    
    ###################
    ## PREPROCESSING ##
    ###################

    # Load matrix of simulated data for training. Consider making this hard-coded in -- path is "/users/j29tien/colocalization_ML/colocalization_matrix_NEW_UNRANKED.tsv"
    matrix_save_path = locus.settings["methods"]["ensemble"]["training_data"] # current matrix contains COLOC, RTC, FINEMAP_STANDARD, FINEMAP_MODIFIED, CAVIARBF, BASELINE, SMARTBASELINE, SMR, SMR_HEIDIADJUSTED, GSMR, TWAS 
    matrix = pd.read_csv(matrix_save_path, sep="\t")

    # remove any submethods based on input and create ScoreContainer numpy array
    methods_used = []
    scores = np.array([])
    if isinstance(locus.scores.coloc_h4,numbers.Number):
	methods_used.append("COLOC")
	scores = np.append(scores,locus.scores.coloc_h4)
    if isinstance(locus.scores.rtc_neg_log_pval,numbers.Number): 
        methods_used.append("RTC")
	scores = np.append(scores,locus.scores.rtc_neg_log_pval)
    if isinstance(locus.scores.finemap_clpp,numbers.Number): 
        methods_used.append("FINEMAP_STANDARD")
        scores = np.append(scores,locus.scores.finemap_clpp)
    if isinstance(locus.scores.finemap_clpp_mod,numbers.Number): 
        methods_used.append("FINEMAP_MODIFIED")
        scores = np.append(scores,locus.scores.finemap_clpp_mod)
    # CAVIARBF not yet implemented
    if isinstance(locus.scores.baseline_neg_log_pval,numbers.Number): 
        methods_used.append("BASELINE")
        scores = np.append(scores,locus.scores.baseline_neg_log_pval)
    if isinstance(locus.scores.smart_baseline_neg_log_pval,numbers.Number): 
        methods_used.append("SMARTBASELINE")
        scores = np.append(scores,locus.scores.smart_baseline_neg_log_pval)
    if isinstance(locus.scores.smr_neg_log_pval,numbers.Number): 
        methods_used.append("SMR")
        scores = np.append(scores,locus.scores.smr_neg_log_pval)
    # SMR_HEIDIADJUSTED not yet implemented
    if isinstance(locus.scores.gsmr_neg_log_pval,numbers.Number): 
        methods_used.append("GSMR")
        scores = np.append(scores,locus.scores.gsmr_neg_log_pval)
    if isinstance(locus.scores.twas_neg_log_pval,numbers.Number): 
        methods_used.append("TWAS")
        scores = np.append(scores,locus.scores.twas_neg_log_pval)
    
    methods_used.append("colocalization_status") # always include answer key
    matrix = matrix[methods_used]

    #print(methods_used)

    # use all data in matrix for training (no need for train test split)
    trainX = matrix.drop('colocalization_status', 1)
    trainY = matrix['colocalization_status']
    
    #convert training and testing data to numpy arrays
    trainX = trainX.to_numpy()  # access the numpy array containing values
    trainY = trainY.to_numpy()  # access the numpy array containing values
    
    # Binarize the output for training
    #trainY = label_binarize(trainY, classes=[0, 1])
    #n_classes = trainY.shape[1]

    ####################
    ## TRAIN ENSEMBLE ##
    ####################
    
    # initiate and fit the random forest classifier -- need to tweak parameters --> n_estimators = 150, max_depth = 7, min_samples_leaf=1
    #randfor = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini', max_depth=7, max_features='auto', max_leaf_nodes=None, min_impurity_decrease=0.0, min_impurity_split=None, min_samples_leaf=1, min_samples_split=2, min_weight_fraction_leaf=0.0, n_estimators=150, n_jobs=None, oob_score=False, random_state=None, verbose=0, warm_start=False)
    
    # for sklearn version < 0.20.0, remove min_impurity_decrease=0.0, set min_impurity_split=0, remove n_jobs=None    
    randfor = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini', max_depth=7, max_features='auto', max_leaf_nodes=None, min_impurity_split=0, min_samples_leaf=1, min_samples_split=2, min_weight_fraction_leaf=0.0, n_estimators=150, oob_score=False, random_state=None, verbose=0, warm_start=False)
    
    randfor.fit(trainX, trainY)
    

    ###############################
    ## COLOCALIZATION PREDICTION ##
    ###############################
    
    #predict using real data
    #print(scores)
    # for sklearn version < 0.20.0, scores.reshape(1,-1) is necessary
    scores = scores.reshape(1,-1)

    predY = randfor.predict(scores)
    #print(predY)
    predProbY = randfor.predict_proba(scores)
    #print("COLOCALIZATION SCORE:")
    #print(predProbY[:,1])
    
    ############
    ## OUTPUT ##
    ############

    # Set ensemble_score to (true) prediction probability 
    # of random forest classifier.
    ensemble_score = predProbY[0,1]

    # Write the names of all the methods / scores that were included (not NoneType) in this file so we can
    # reference it later
    methods_used.remove('colocalization_status')
    with open("{0}/ensemble_methods_used.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as w:
        w.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, locus.data.shape[0], locus.gwas_suffix))
	for method in methods_used:
            w.write("%s," % method)       
	w.write("\n") 

    with open("{0}/{1}_ensemble_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, locus.data.shape[0], locus.gwas_suffix, ensemble_score))

    return ensemble_score
