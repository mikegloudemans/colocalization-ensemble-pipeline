#!/usr/bin/python
# Author: Mike Gloudemans
# 
# Dispatch colocalization analyses
# for various GWAS and eQTL experiments,
# parameter settings, and colocalization
# methods.
#

import subprocess
import pandas as pd
import sys
import operator
import time
import argparse
import os

# Definitely need the libraries below here.
# Above -- not so sure.

# Built-ins
from shutil import copyfile

# Custom
import config
import preprocess
import tabix_snps
import TestLocus

def main():

        # TODO: Make results directory, under which all output for this run
        # will be stored.
        # It must be timestamped in an intuitive way, unique to this particular run.
        # We will pass this to other classes via the Experiment object.
        # Probably should be a separate function.
	base_output_dir = "" # TODO
        # TODO: Make it using a shell command

        # Hard-coded for now; will add an argparse function to do this later.
        config_file = "/users/mgloud/projects/brain_gwas/data/config/sample.config"

        # Read config file
        settings = config.load_config(config_file)
        # Copy config file contents to the output directory so we know what
        # parameters were used for this run.
        copyfile(config_file, "{0}/{1}/settings_used.config".format(base_output_dir, "/".join(config_file.split("/")[:-1])))
       
        # Index all SNP files as needed.
        # For now we're not going to index GWAS SNPs because they're not
        # big enough to need this.
        tabix_snps(settings)

        # For each GWAS experiment:
        for gwas_file in settings.gwas_files:

            # Call select_test_snps() (in another module) to
            # get a list of which SNPs we should test in this GWAS.
            snp_list = preprocess.select_test_snps(gwas_file)

            # For each eQTL experiment:
            for eqtl_file in settings.eqtl_files:

                # For each GWAS SNP selected above...
                for snp in snp_list:

                    # Load relevant GWAS and eQTL data using tabix.
                    # Merge the two tables using a helper class with pandas,
                    # similar to what we were doing before.
                    # (Separate module -- get the resulting pandas table.)
                    data = preprocess.load_summary_statistics(snp, gwas_file, eqtl_file)

                    # Create a TestLocus object using merged GWAS and eQTL,
                    # any important metadata about the experiment such as the directory,
                    # and the Config object.
                    task = TestLocus(data, settings, base_output_dir)
                    task.run()


if __name__ == "__main__":
	main()
        
