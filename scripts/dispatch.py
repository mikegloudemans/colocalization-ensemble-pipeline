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
from multiprocessing import Pool
import traceback
import gzip
import os

# Custom libraries
import config
import preprocess
from TestLocus import TestLocus

def main():

    # Read config file
    config_file = sys.argv[1]
    settings = config.load_config(config_file)

    max_cores = int(sys.argv[2])

    # Make timestamped results directory, under which all output for this run will be stored.
    # Note: sometimes this might conflict with another run of the script. If so, keep trying until
    # a directory name is free
    while True:
        now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S.%f')
        if not os.path.isdir("/users/mgloud/projects/brain_gwas/tmp/{0}".format(now)):  
            try:
                subprocess.check_call("mkdir /users/mgloud/projects/brain_gwas/tmp/{0}".format(now), shell=True)
                # If another script beats us to it, this we'll fail so we'll try again then
            except:
                continue
            break


    if "out_dir_group" in settings:
        base_output_dir = "/users/mgloud/projects/brain_gwas/output/{0}/{1}".format(settings["out_dir_group"], now)
    else:
        base_output_dir = "/users/mgloud/projects/brain_gwas/output/{0}".format(now)
    base_output_dir = base_output_dir + "_" + config_file.split("/")[-1].split(".")[0]
    base_tmp_dir = "/users/mgloud/projects/brain_gwas/tmp/{0}".format(now)

    # Save config file and current Git log for reproducibility.
    save_state(config_file, base_output_dir)

    gwas_files = sorted([f for f in settings["gwas_experiments"]])
    eqtl_files = sorted([f for f in settings["eqtl_experiments"]])

    # For each GWAS experiment:
    for gwas_file in gwas_files:

        # Write header of output file for FINEMAP

        gwas_suffix = gwas_file.split("/")[-1].replace(".", "_") 
        if "finemap" in settings["methods"]:
            with open("{0}/{1}_finemap_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tclpp\t-log_gwas_pval\t-log_eqtl_pval\tbase_gwas_file\tclpp_mod\n")

        # Write COLOC results to the desired file.
        if "coloc" in settings["methods"]:
            with open("{0}/{1}_coloc_h4pp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tclpp_h4\tbase_gwas_file\n")

        if "ecaviar" in settings["methods"]:
            with open("{0}/{1}_ecaviar_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tclpp\n")

        if "rtc" in settings["methods"]:
            with open("{0}/{1}_rtc_score_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\ttrait\tfeature\trtc_score\tbase_gwas_file\n")
 
        if "caviarbf" in settings["methods"]:
            with open("{0}/{1}_caviarbf_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tclpp\n")

        if "twas" in settings["methods"]:
            with open("{0}/{1}_twas_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\ttwas_log_pval\ttwas_perm_log_pval\n")

        if "metaxcan" in settings["methods"]:
            with open("{0}/{1}_metaxcan_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\ttwas_log_pval\n")

        if "baseline" in settings["methods"]:
            with open("{0}/{1}_baseline_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tbase_gwas_file\tbaseline_pval\tbaseline_pval2\tbaseline_pval3\tbaseline_pval4\tbaseline_pval5\n")



        # Get list of traits measured in this GWAS
        traits = set([])

        # Are there multiple traits in this GWAS?
        if "traits" not in settings["gwas_experiments"][gwas_file]:
            with gzip.open(gwas_file) as f:
                header = f.readline().strip().split("\t")
                if "trait" in header:
                    trait_index = header.index("trait")
                    for line in f:
                        traits.add(line.strip().split("\t")[trait_index])
                else:
                    traits.add(gwas_file.split("/")[-1])
            traits = list(traits)

        # Subset down to traits of interest, if specified
        if "traits" in settings["gwas_experiments"][gwas_file]:
            traits = settings["gwas_experiments"][gwas_file]["traits"]

        assert len(traits) != 0

        for trait in traits:

            print "Testing", trait

            gwas_snp_list = []
            # Get a list of which SNPs we should test in this GWAS.

            # Selection basis options:
            #   - gwas: SNPs significant in GWAS will be tested (at specified threshold)
            #   - eqtl: SNPs significant in eQTLs will be tested (at specified threshold)
            #   - both: SNPs significant in GWAS or eQTL will be tested
            # To run for the entire genome, specify "eQTL" and set the pvalue cutoff to 1.

            if settings["selection_basis"] in ["gwas", "both"]:
                gwas_snp_list.extend(preprocess.select_test_snps_by_gwas(gwas_file, settings['selection_thresholds']["gwas"], trait))

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
                    pool.apply_async(analyze_snp_wrapper, args=(gwas_file, eqtl_file, snp[0], settings, base_output_dir, base_tmp_dir, trait), kwds=dict(restrict_gene=snp[1]))
                pool.close()
                pool.join()
     
                # Clean up after ourselves
                subprocess.call("rm -r {0} 2> /dev/null".format(base_tmp_dir), shell=True)

                # Run GWAS SNPs separately just in case there happen to be any overlaps,
                # which could lead to a race.
                pool = Pool(max_cores)
                for i in xrange(0, len(gwas_snp_list)):
                    snp = gwas_snp_list[i]
                    pool.apply_async(analyze_snp_wrapper, args=(gwas_file, eqtl_file, snp[0], settings, base_output_dir, base_tmp_dir, trait), kwds=dict(restrict_gene=snp[1]))
                pool.close()
                pool.join()

                # Clean up after ourselves
                subprocess.call("rm -r {0} 2> /dev/null".format(base_tmp_dir), shell=True)

                # Make SplicePlots if appropriate
                if "splice_plots" in settings and eqtl_file in settings["splice_plots"]:
                    results_file = "{0}/{1}_finemap_clpp_status.txt".format(base_output_dir, gwas_suffix)
                    splice_plot(results_file, eqtl_file, settings)

           
    # Create full genome-wide plot of results (currently just for CLPP - TODO fix)
    # TODO: Move this to a separate function
    #for gwas_file in gwas_files:
    #    gwas_suffix = gwas_file.split("/")[-1].replace(".", "_")
    #
    #    subprocess.check_call("mkdir -p {0}/manhattan".format(base_output_dir), shell=True)
    #    subprocess.check_call("Rscript /users/mgloud/projects/brain_gwas/scripts/full_genome_plot.R {0}/{1}_finemap_clpp_status.txt {0}/manhattan".format(base_output_dir, gwas_suffix), shell=True)


# If we're running in parallel and a thread fails, catch the exception
# and log it to a file so we can fix it later.
def analyze_snp_wrapper(gwas_file, eqtl_file, snp, settings, base_output_dir, base_tmp_dir, trait, restrict_gene=-1):
    try:
        analyze_snp(gwas_file, eqtl_file, snp, settings, base_output_dir, base_tmp_dir, trait, restrict_gene)
    except Exception as e:
        print "caught exception"
        traceback.print_exc(file=sys.stdout)
        error = str(e)
        error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
        with open("{0}/ERROR_variants.txt".format(base_output_dir),"a") as a:
            a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, restrict_gene, trait, error))
        raise Exception("Failed coloc run.")

def analyze_snp(gwas_file, eqtl_file, snp, settings, base_output_dir, base_tmp_dir, trait, restrict_gene=-1):

    # Load relevant GWAS and eQTL data.
    gwas_data = preprocess.get_gwas_data(gwas_file, snp, settings, trait) # Get GWAS data
    eqtl_data = preprocess.get_eqtl_data(eqtl_file, snp, settings) # Get eQTL data

    # Skip it if no GWAS variants at this locus
    if isinstance(gwas_data, basestring):
        # Write skipped variants to a file, for later reference.
        with open("{0}/skipped_variants.txt".format(base_output_dir),"a") as a:
            a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, "-1", gwas_data, trait))
        return

    # Skip it if this entire locus has no genes
    if isinstance(eqtl_data, basestring):
        # Write skipped variants to a file, for later reference.
        with open("{0}/skipped_variants.txt".format(base_output_dir),"a") as a:
            a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, "-1", eqtl_data, trait))
        return

    # Temporary mod for splice eQTLs. May be best to specify a "feature" ID in the future
    if 'feature' in eqtl_data:
        eqtl_data['gene'] = eqtl_data['feature']
    if 'Gene' in eqtl_data:
        eqtl_data['gene'] = eqtl_data['Gene']
    # Don't want to have colons in our filenames later
    eqtl_data['gene'] = eqtl_data['gene'].str.replace(':', '.')
    # Don't want our temporary filename to be too long
    def trim_trait(x):
        if len(x) > 100:
            return x[:100] + "_trimmed"
        else:
            return x
    eqtl_data['gene'] = eqtl_data['gene'].apply(trim_trait)

    # Get all genes whose eQTLs we're testing at this locus
    if isinstance(restrict_gene, basestring) and len(restrict_gene) > 100:
        restrict_gene = restrict_gene[:100] + "_trimmed"
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
                a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, gene, combined, trait))
            continue

        # Create a TestLocus object using merged GWAS and eQTL,
        # any important metadata about the experiment such as the directory,
        # and the Config object.
        task = TestLocus(combined, settings, base_output_dir, base_tmp_dir, gene, snp, gwas_file, eqtl_file, trait)
        task.run()

# At start of run, save settings so we'll know what they were when we ran it.
def save_state(config_file, base_output_dir):

    # Copy config file to output for later reference
    subprocess.check_call("mkdir -p {0}".format(base_output_dir), shell=True)
    copyfile(config_file, "{0}/settings_used.config".format(base_output_dir))

    # For reproducibility, store the current state of the project in Git
    subprocess.check_call('git --git-dir /users/mgloud/projects/brain_gwas/.git log >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    subprocess.check_call('git --git-dir /users/mgloud/projects/brain_gwas/.git diff >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    subprocess.check_call('git --git-dir /users/mgloud/projects/brain_gwas/.git branch >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    subprocess.check_call('git --git-dir /users/mgloud/projects/brain_gwas/.git status >> {0}/git_status.txt'.format(base_output_dir), shell=True)
 
# Use SplicePlot to generate splice plots
def splice_plot(results_file, eqtl_file, settings):
    map_file = settings["splice_plots"][eqtl_file]["map_file"]
    vcf_file = settings["splice_plots"][eqtl_file]["vcf_file"]
    sp.generate_splice_plots(results_file, eqtl_file, vcf_file, map_file, threshold = 0.01)

if __name__ == "__main__":
    main()
 
