#!/usr/bin/python
# Author: Mike Gloudemans
# 
# Dispatch colocalization analyses
# for various GWAS and eQTL experiments,
# parameter settings, and colocalization
# methods.
#

max_cores = 12

# Built-in libraries
import sys
from shutil import copyfile
import datetime
import subprocess
import math
from multiprocessing import Pool
import traceback

# Custom libraries
import config
import preprocess
from TestLocus import TestLocus

def main():

    # Read config file
    config_file = sys.argv[1]
    settings = config.load_config(config_file)
 
    # Make timestamped results directory, under which all output for this run will be stored.
    now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    base_output_dir = "/users/mgloud/projects/brain_gwas/output/{0}".format(now)
    base_output_dir = base_output_dir + "_" + config_file.split("/")[-1].split(".")[0]
    base_tmp_dir = "/users/mgloud/projects/brain_gwas/tmp/{0}".format(now)

    # Save config file and current Git log for reproducibility.
    save_state(config_file, base_output_dir)

    gwas_files = [f for f in settings["gwas_experiments"]]
    eqtl_files = [f for f in settings["eqtl_experiments"]]

    # For each GWAS experiment:
    for gwas_file in gwas_files:

        # Write header of output file for FINEMAP
        gwas_suffix = gwas_file.split("/")[-1].replace(".", "_") 
        with open("{0}/{1}_finemap_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
            w.write("ref_snp\teqtl_file\tgwas_file\tfeature\tn_snps\tclpp\t-log_gwas_pval\t-log_eqtl_pval\n")

        gwas_snp_list = []
        # Get a list of which SNPs we should test in this GWAS.

        # Selection basis options:
        #   - gwas: SNPs significant in GWAS will be tested (at specified threshold)
        #   - eqtl: SNPs significant in eQTLs will be tested (at specified threshold)
        #   - both: SNPs significant in GWAS or eQTL will be tested
        # To run for the entire genome, specify "eQTL" and set the pvalue cutoff to 1.

        if settings["selection_basis"] in ["gwas", "both"]:
            gwas_snp_list.extend(preprocess.select_test_snps_by_gwas(gwas_file, settings['gwas_threshold']))

        if settings["selection_basis"] == "snps_from_list":
            gwas_snp_list.extend(preprocess.select_snps_from_list(settings["snp_list_file"]))

        # For each eQTL experiment:
        for eqtl_file in eqtl_files:

            eqtl_snp_list = []
            if settings["selection_basis"] in ["eqtl", "both"]:
                # If a "selection subset" is specified for the eQTL experiment, then genes will
                # only be tested if they are in this subset.
                if "selection_subset" in settings['eqtl_experiments'][eqtl_file]:
                    eqtl_subset = settings['eqtl_experiments'][eqtl_file]["selection_subset"]
                else:
                    eqtl_subset = -1
                eqtl_snp_list.extend(preprocess.select_test_snps_by_eqtl(eqtl_file, settings, eqtl_subset))

            snp_list = eqtl_snp_list + gwas_snp_list
            print("Testing {2} SNPs ({0} GWAS hits and {1} eQTL hits).".format(len(gwas_snp_list), len(eqtl_snp_list), len(snp_list)))

            # Run key SNPs in parallel
            pool = Pool(max_cores)
            for i in xrange(0, len(eqtl_snp_list)):
                snp = eqtl_snp_list[i]
                pool.apply_async(analyze_snp_wrapper, args=(gwas_file, eqtl_file, snp[0], settings, base_output_dir, base_tmp_dir), kwds=dict(restrict_gene=snp[1]))
            pool.close()
            pool.join()
 
            # Clean up after ourselves
            subprocess.call("rm -r {0}".format(base_tmp_dir), shell=True)

            # Run GWAS SNPs separately just in case there happen to be any overlaps,
            # which could lead to a race.
            pool = Pool(max_cores)
            for i in xrange(0, len(gwas_snp_list)):
                snp = gwas_snp_list[i]
                pool.apply_async(analyze_snp_wrapper, args=(gwas_file, eqtl_file, snp[0], settings, base_output_dir, base_tmp_dir), kwds=dict(restrict_gene=snp[1]))
            pool.close()
            pool.join()

            # Clean up after ourselves
            subprocess.call("rm -r {0}".format(base_tmp_dir), shell=True)

            # Make SplicePlots if appropriate
            if "splice_plots" in settings and eqtl_file in settings["splice_plots"]:
                results_file = "{0}/{1}_finemap_clpp_status.txt".format(base_output_dir, gwas_suffix)
                splice_plot(results_file, eqtl_file, settings)

           
    # Create full genome-wide plot of results (currently just for CLPP - TODO fix)
    # TODO: Move this to a separate function
    for gwas_file in gwas_files:
        gwas_suffix = gwas_file.split("/")[-1].replace(".", "_")

        subprocess.check_call("mkdir -p {0}/manhattan".format(base_output_dir), shell=True)
        subprocess.check_call("Rscript /users/mgloud/projects/brain_gwas/scripts/full_genome_plot.R {0}/{1}_finemap_clpp_status.txt {0}/manhattan".format(base_output_dir, gwas_suffix), shell=True)


# If we're running in parallel and a thread fails, catch the exception
# and log it to a file so we can fix it later.
def analyze_snp_wrapper(gwas_file, eqtl_file, snp, settings, base_output_dir, base_tmp_dir, restrict_gene=-1):
    try:
        analyze_snp(gwas_file, eqtl_file, snp, settings, base_output_dir, base_tmp_dir, restrict_gene)
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        with open("{0}/ERROR_variants.txt".format(base_output_dir),"a") as a:
            a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, restrict_gene, str(e)))

def analyze_snp(gwas_file, eqtl_file, snp, settings, base_output_dir, base_tmp_dir, restrict_gene=-1):

    # Load relevant GWAS and eQTL data.
    gwas_data = preprocess.get_gwas_data(gwas_file, snp, settings) # Get GWAS data
    eqtl_data = preprocess.get_eqtl_data(eqtl_file, snp, settings) # Get eQTL data

    # Skip it if this entire locus has no genes
    if isinstance(eqtl_data, basestring):
        # Write skipped variants to a file, for later reference.
        with open("{0}/skipped_variants.txt".format(base_output_dir),"a") as a:
            a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, "-1", eqtl_data))
        return

    # Temporary mod for splice eQTLs. May be best to specify a "feature" ID in the future
    if 'feature' in eqtl_data:
        eqtl_data['gene'] = eqtl_data['feature']
    # Don't want to have colons in our filenames later
    eqtl_data['gene'] = eqtl_data['gene'].str.replace(':', '.')

    # Get all genes whose eQTLs we're testing at this locus
    if restrict_gene == -1:
        genes = set(eqtl_data['gene'])
    else:
        restrict_gene_mod = restrict_gene.replace(":", ".") 
        genes = [restrict_gene_mod]

    # Loop through all genes now
    for gene in genes:
        # NOTE: It might be easier to just do this step once outside of this loop,
        # and then filter down to the gene of interest. Consider modifying.
        allow_insignificant_gwas = restrict_gene != -1
        combined = preprocess.combine_summary_statistics(gwas_data, eqtl_data, gene, snp, settings, unsafe=True, allow_insignificant_gwas=allow_insignificant_gwas)

        # Skip it if this site is untestable.
        if isinstance(combined, basestring):
            # Write skipped variants to a file, for later reference.
            with open("{0}/skipped_variants.txt".format(base_output_dir),"a") as a:
                a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, gene, combined))
            continue

        # Create a TestLocus object using merged GWAS and eQTL,
        # any important metadata about the experiment such as the directory,
        # and the Config object.
        task = TestLocus(combined, settings, base_output_dir, base_tmp_dir, gene, snp, gwas_file, eqtl_file)
        task.run()

# At start of run, save settings so we'll know what they were when we ran it.
def save_state(config_file, base_output_dir):

    # Copy config file to output for later reference
    subprocess.check_call("mkdir -p {0}".format(base_output_dir), shell=True)
    copyfile(config_file, "{0}/settings_used.config".format(base_output_dir))

    # For reproducibility, store the current state of the project in Git
    subprocess.check_call('git log >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    subprocess.check_call('git diff >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    subprocess.check_call('git branch >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    subprocess.check_call('git status >> {0}/git_status.txt'.format(base_output_dir), shell=True)
 
# Use SplicePlot to generate splice plots
def splice_plot(results_file, eqtl_file, settings):
    map_file = settings["splice_plots"][eqtl_file]["map_file"]
    vcf_file = settings["splice_plots"][eqtl_file]["vcf_file"]
    sp.generate_splice_plots(results_file, eqtl_file, vcf_file, map_file, threshold = 0.01)

if __name__ == "__main__":
	main()
 
