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
import math

# Custom libraries
import config
import preprocess
import tabix_snps
from TestLocus import TestLocus

# TODO: Parallelize this again. As it is right now, it's fairly slow because so many sites/tissues to test.

# TODO: Make TMP files be time-stamped so we never have to worry about conflicts again

def main():

    config_file = sys.argv[1]

    # Make timestamped results directory, under which all output for this run will be stored.
    now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    base_output_dir = "/users/mgloud/projects/brain_gwas/output/{0}".format(now)


    # Read config file
    settings = config.load_config(config_file)
    if "output_suffix" in settings:
        base_output_dir = base_output_dir + "_" + settings["output_suffix"]
    else:
        base_output_dir = base_output_dir + "_" + config_file.split("/")[-1].split(".")[0]

    # Copy config file to output for later reference
    subprocess.check_call("mkdir -p {0}".format(base_output_dir), shell=True)
    copyfile(config_file, "{0}/settings_used.config".format(base_output_dir))

    # For reproducibility, store the current state of the project in Git
    subprocess.check_call('git log >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    subprocess.check_call('git diff >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    subprocess.check_call('git branch >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    subprocess.check_call('git status >> {0}/git_status.txt'.format(base_output_dir), shell=True)

    # TODO: Index all SNP files as needed, for both GWAS and eQTLs.
    # Up to this point we've just been doing it manually.
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
                gwas_data = preprocess.get_gwas_data(gwas_file, snp, settings) # Get GWAS data
                eqtl_data = preprocess.get_eqtl_data(eqtl_file, snp, settings) # Get eQTL data

                # Skip it if this entire locus has no genes
                if isinstance(eqtl_data, basestring):
                    # Write skipped variants to a file, for later reference.
                    with open("{0}/skipped_variants.txt".format(base_output_dir),"a") as a:
                        a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, "-1", eqtl_data))
                    continue

                # Temporary mod for splice eQTLs. May be best to specify a "feature" ID in the future
                if 'feature' in eqtl_data:
                    eqtl_data['gene'] = eqtl_data['feature']

                # Get all genes whose eQTLs we're testing at this locus
                genes = set(eqtl_data['gene'])
               
                # Loop through all genes now
                for gene in genes:
                    # NOTE: It might be easier to just do this step once outside of this loop,
                    # and then filter down to the gene of interest. Consider modifying.
                    combined = preprocess.combine_summary_statistics(gwas_data, eqtl_data, gene, snp, settings, unsafe=True)

                    # Skip it if this site is untestable.
                    if isinstance(combined, basestring):
                        # Write skipped variants to a file, for later reference.
                        with open("{0}/skipped_variants.txt".format(base_output_dir),"a") as a:
                            a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, gene, combined))
                        continue

                    # Create a TestLocus object using merged GWAS and eQTL,
                    # any important metadata about the experiment such as the directory,
                    # and the Config object.
                    task = TestLocus(combined, settings, base_output_dir, gene, snp, gwas_file, eqtl_file)
                    task.run()

    # Create full genome-wide plot of results (currently just for CLPP - TODO fix)
    # TODO: Move this to a separate function
    for gwas_file in gwas_files:
        gwas_suffix = gwas_file.split("/")[-1].split(".")[0].replace(".", "_")

        subprocess.check_call("mkdir -p {0}/manhattan", shell=True)
        subprocess.check_call("Rscript /users/mgloud/projects/brain_gwas/scripts/full_genome_plot.R {0}/{1}_finemap_clpp_status.txt {0}/manhattan".format(base_output_dir, gwas_suffix), shell=True)

    # Clean up after ourselves

    # Don't do this. This is dangerous because if multiple instances of the pipeline are running
    # at the same time, one can erase the other's data.
    #subprocess.call("rm -r /users/mgloud/projects/brain_gwas/tmp/*", shell=True)

if __name__ == "__main__":
	main()
        
