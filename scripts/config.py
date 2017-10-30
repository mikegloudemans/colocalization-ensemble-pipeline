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

    return config
