#!/usr/bin/python
#
# Author: Mike Gloudemans
# 
# Tabix lists of SNPs for easy access.
# 

import subprocess
import os

# Function tabix_file
# 
# Arguments: Config object
# Action: tabix all files in Config.
# Returns: None
def tabix_all(config):

    # Call tabix_file for each file in config.
    for file in [f for f in config["eqtl_experiments"]]:
        tabix_file(file)


# Arguments: name of file, whether to reindex if file exists already
# Action: tabix the file
# Returns: None
def tabix_file(filename, redo=False):

    # Check whether file has already been tabixed. If so,
    # quit.
    if os.path.isfile("{0}.tbi".format(filename)):
        return

    # Otherwise, tabix the file.
    print("Indexing files with tabix. This may take several hours, but only is necessary \
            on the first run with a given file.")

    # Run tabix
    # TODO: Test to see if this code words. Right now untested.

    #subprocess.check_call("gunzip {0}".format(filename), shell=True)
    #base = filename.split(".gz")[0]
    #subprocess.check_call('''cat <(echo "chr snp_pos ref     alt     genome  gene    beta    t-stat  pvalue") <(tail -n +2 {0} | sed 's/_/\\t/g' | sort -k1,1n -k2,2n) | bgzip > {0}.gz'''.format(base) shell=True)
    #subprocess("tabix $base.gz -b 2 -e 2 -s 1 -S 1".format(base), shell=True)
    # TODO: Unsafe: consider unzip to a temp folder so we don't accidentally erase the original
    # file if the user wanted it.
    #subprocess.check_call("rm {0}".format(base), shell=True)
