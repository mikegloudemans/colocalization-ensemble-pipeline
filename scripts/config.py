# Author: Mike Gloudemans
# 
# Load and parse config file.
#

import json

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


    return config
