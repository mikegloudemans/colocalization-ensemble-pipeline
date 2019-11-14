# Author: Mike Gloudemans
# 
# Load and parse config file.
#

import json 
import glob
import copy

debug = False

# Function load_config
#
# Input: Filename
# Returns: dictionary containing config settings.
def load_config(filename):

    with open(filename) as data_file:    
            config = json.load(data_file)

    if "eqtl_threshold" not in config:
        config["eqtl_threshold"] = 1e-7

    # Default to effect size format if no format specified
    for file in config['eqtl_experiments']:
        if 'eqtl_format' not in config['eqtl_experiments'][file]:
            config['eqtl_experiments'][file]['eqtl_format'] = "effect_size"

    for file in config['gwas_experiments']:
        if 'gwas_format' not in config['gwas_experiments'][file]:
            config['gwas_experiments'][file]['gwas_format'] = "effect_size"

    if "plot_all" in config and config["plot_all"]=="True":
        config["plot_all"] = True
    else:
        config["plot_all"] = False
    
    if "plot_none" in config and config["plot_none"]=="True":
        config["plot_none"] = True
    else:
        config["plot_none"] = False

    if "plot_only" in config and config["plot_only"]=="True":
        config["plot_only"] = True
    else:
        config["plot_only"] = False

    # Default window of analysis = 500000 bp on either side of SNP
    if "window" not in config:
        config["window"] = 500000

    if debug:
        config["debug"] = True
    if "debug" in config and config["debug"] == "True":
        config["debug"] = True
    else:
        config["debug"] = False

    config = expand_glob_eqtl_files(config)

    for ee in config['eqtl_experiments']:
        if "selection_subset" in config['eqtl_experiments'][ee]:
            config['eqtl_experiments'][ee]["selection_subset"] = get_selection_subset(config['eqtl_experiments'][ee]["selection_subset"])

    return config

def get_selection_subset(filename):
    selection_subset = set([])
    with open(filename) as f:
        for line in f:
            selection_subset.add(line.strip())
    return selection_subset

def expand_glob_eqtl_files(config):

    new_config = copy.deepcopy(config)

    for ee in config['eqtl_experiments']:
        all_experiments = glob.glob(ee)
        del new_config['eqtl_experiments'][ee]

        for ae in all_experiments:
            new_config['eqtl_experiments'][ae] = copy.deepcopy(config['eqtl_experiments'][ee])
    return new_config
