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

    return config

# For testing purposes only:
if __name__ == "__main__":
    c = load_config("/users/mgloud/projects/brain_gwas/data/config/sample.config")
    print c
