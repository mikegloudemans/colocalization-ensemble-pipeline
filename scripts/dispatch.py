#!/usr/bin/python
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

# Custom libraries
import config
import preprocess
import tabix_snps
from TestLocus import TestLocus

def main():

    # Hard-coded for now; will add an argparse function to do this later.
    config_file = "/users/mgloud/projects/brain_gwas/data/config/sample.config"

    # Make timestamped results directory, under which all output for this run will be stored.
    now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    base_output_dir = "/users/mgloud/projects/brain_gwas/output/{0}".format(now)
    subprocess.check_call("mkdir -p {0}".format(base_output_dir), shell=True)

    # Read config file
    settings = config.load_config(config_file)
    # Copy config file to output for later reference
    copyfile(config_file, "{0}/settings_used.config".format(base_output_dir))
   
    # Index all SNP files as needed.
    # For now we're not going to index GWAS SNPs because they're not
    # big enough to need this.
    tabix_snps.tabix_all(settings)

    gwas_files = [f for f in settings["gwas_experiments"]]
    eqtl_files = [f for f in settings["eqtl_experiments"]]

    # For each GWAS experiment:
    for gwas_file in gwas_files:

        # Get a list of which SNPs we should test in this GWAS.
        snp_list = preprocess.select_test_snps(gwas_file, settings["gwas_threshold"])

        # For each eQTL experiment:
        for eqtl_file in eqtl_files:

            # For each GWAS SNP selected above...
            for snp in snp_list:

                # Load relevant GWAS and eQTL data.
                gwas_data = preprocess.get_gwas_data(gwas_file, snp) # Get GWAS data
                eqtl_data = preprocess.get_eqtl_data(eqtl_file, snp) # Get eQTL data

                # Get full all genes whose eQTLs we're testing at this locus
                genes = set(eqtl_data['gene'])
               
                # Loop through all genes now
                for gene in genes:

                    combined = preprocess.combine_summary_statistics(gwas_data, eqtl_data, gene, snp)
                    # Skip it if this site is untestable.
                    if combined is None:
                        continue

                    # Create a TestLocus object using merged GWAS and eQTL,
                    # any important metadata about the experiment such as the directory,
                    # and the Config object.
                    task = TestLocus(combined, settings, base_output_dir, gene, snp, gwas_file, eqtl_file)
                    task.run()


    # TODO: Post-analysis of results.

if __name__ == "__main__":
	main()
        
