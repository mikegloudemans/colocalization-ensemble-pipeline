#!/usr/bin/python3
# Author: Mike Gloudemans
#
# Dispatch colocalization analyses
# for various GWAS and eQTL experiments,
# parameter settings, and colocalization
# methods.
#

# Built-in libraries
import sys
from shutil import copyfile
import datetime
import subprocess
import math
from multiprocessing import Pool
import traceback
import gzip
import os
import time 
from progress.bar import Bar
import json
import re

# Custom libraries
import config
import preprocess
from TestLocus import TestLocus
import dispatch

# do some things before we start the pipeline 
# Read config file
config_file = "/Users/nicolerg/Downloads/page_config.json"
# Check if absolute path or relative path
if not config_file.startswith("/"):
    config_file = os.getcwd() + "/" + config_file

# Change to directory of script
os.chdir(os.path.abspath(os.path.dirname(sys.argv[0])))

#settings = config.load_config(config_file)
with open(config_file) as data_file:    
    settings = json.load(data_file)

# verify that all GWAS and eQTL files exist; if not, abort with an error message.
for f in settings["gwas_experiments"]:
    if not os.path.exists(f):
        raise Exception("Error: requested GWAS file {0} does not exist.".format(f))
for f in settings["eqtl_experiments"]:
    if not os.path.exists(f):
        raise Exception("Error: requested eQTL file {0} does not exist.".format(f))

# check that out_dir and tmp_dir are specified and exist
if not "out_dir" in settings:
    raise Exception("Error: 'out_dir' not specified in config file.")
else: 
    out_dir = settings["out_dir"]
if not os.path.exists(out_dir):	
    raise Exception("Error: specified 'out_dir' DNE: {0}".format(out_dir))

if not "tmp_dir" in settings:
    raise Exception("Error: 'tmp_dir' not specified in config file.")
else: 
    tmp_dir = settings["tmp_dir"]
if not os.path.exists(tmp_dir):	
    raise Exception("Error: specified 'tmp_dir' DNE: {0}".format(tmp_dir))

# Make timestamped results directory, under which all output for this run will be stored.
# Note: sometimes this might conflict with another run of the script. If so, keep trying until
# a directory name is free
while True:
    now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S.%f')
    if not os.path.isdir("{0}/{1}".format(tmp_dir, now)):
        try:
            subprocess.check_call("mkdir {0}/{1}".format(tmp_dir, now), shell=True)
            # If another script beats us to it, this we'll fail so we'll try again then
        except:
            continue
        break

if "out_dir_group" in settings:
    base_output_dir = "{0}/{1}/{2}".format(out_dir, settings["out_dir_group"], now)
else:
    base_output_dir = "{0}/{1}".format(out_dir, now)
base_output_dir = base_output_dir + "_" + config_file.split("/")[-1].split(".")[0]
base_tmp_dir = "{0}/{1}".format(tmp_dir, now)

# Save config file and current Git log for reproducibility.
dispatch.save_state(config_file, base_output_dir)

# trait names, not file names
gwas_list = sorted([f for f in settings["gwas_experiments"]])
eqtl_list = sorted([f for f in settings["eqtl_experiments"]])
methods = sorted([m for m in settings["methods"]])

# add "trait" attribute to gwas 
for gwas in gwas_list:

	# get file name
	gwas_file = settings["gwas_experiments"][gwas]["file"]

	# Get list of traits measured in this GWAS
    traits = set([])
    # Are there multiple traits in this GWAS?
    # if not, assign "traits" == gwas 
    if "traits" not in settings["gwas_experiments"][gwas]:
        with gzip.open(gwas_file) as f:
            header = f.readline().strip().split("\t")
            if "trait" in header:
                trait_index = header.index("trait")
                for line in f:
                    traits.add(line.strip().split("\t")[trait_index])
                # add traits to settings
                traits = list(traits)
                settings["gwas_experiments"][gwas]["traits"] = traits
            else:
            	settings["gwas_experiments"][gwas]["traits"] = gwas

traits = list()

for gwas in settings["gwas_experiments"]:
	for trait in settings["gwas_experiments"][gwas]:
		traits.add(trait)



# what do I want in the end?
# one file per gwas-eqtl file pair per method 
rule all:
    input:
        expand("{eqtl}.{gwas}.{method}.results.tsv", eqtl=eqtl_list, gwas=gwas_list, method=methods)



	





