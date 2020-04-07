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
import time 
from progress.bar import Bar
import pandas as pd
import logging

# Custom libraries
import config
import preprocess
from TestLocus import TestLocus
import SNP

def main():

    # add time stamps to log file 
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    # Read config file
    logging.info('Load settings...')
    config_file = sys.argv[1]
    # Check if absolute path or relative path
    if not config_file.startswith("/"):
       config_file = os.getcwd() + "/" + config_file

    # Change to directory of script
    os.chdir(os.path.abspath(os.path.dirname(sys.argv[0])))

    settings = config.load_config(config_file)
	
    logging.info('Verify that input files exist...')
    if not "overlap_loci" in settings["selection_basis"]:
        
	# Verify that all GWAS and eQTL files exist; if not, abort with an error message.
        for f in settings["gwas_experiments"]:
            if not os.path.exists(f):
                raise Exception("Error: requested GWAS file {0} does not exist.".format(f))
        for f in settings["eqtl_experiments"]:
            if not os.path.exists(f):
                raise Exception("Error: requested eQTL file {0} does not exist.".format(f))
		
    else:
	
	# check for existence of files present in "overlap_loci" file
	gwas_files = set([])
	eqtl_files = set([])
	
	# check required columns and existence of files 
	with open(settings["selection_basis"]["overlap_loci"], 'r') as f:
	    
	    header = f.readline().strip().split("\t")
	    
	    check_columns = ['chr' in header,
			     'snp_pos' in header,
			     'gwas_pvalue' in header,
			     'eqtl_pvalue' in header,
			     'feature' in header,
			     'gwas_file' in header,
			     'eqtl_file' in header]
	
	    if not all(check_columns):
		raise Exception("Error: A required column in {} is missing. Please see 'README.txt.'.".format(settings["selection_basis"]["overlap_loci"]))
	
	    gwas_index = header.index("gwas_file")
	    eqtl_index = header.index("eqtl_file")
	    for line in f:
		gwas_files.add(line.strip().split("\t")[gwas_index])
		eqtl_files.add(line.strip().split("\t")[eqtl_index])
	
	for f in gwas_files:
	    if not os.path.exists(f):
	        raise Exception("Error: requested GWAS file {0} does not exist.".format(f))
	for f in eqtl_files:
	    if not os.path.exists(f):
		raise Exception("Error: requested eQTL file {0} does not exist.".format(f))
			
    max_cores = int(sys.argv[2]) # this is actually threads, not cores 

    logging.info('Set up directories...')
    if "out_dir" in settings:
        out_dir = settings["out_dir"]
    else:
        raise Exception("Error: 'out_dir' not specified.")

    if "tmp_dir" in settings:
        tmp_dir = settings["tmp_dir"]
    else:
        raise Exception("Error: 'tmp_dir' not specified.")

    # Make timestamped results directory, under which all output for this run will be stored.
    # Note: sometimes this might conflict with another run of the script. If so, keep trying until
    # a directory name is free
    while True:
        now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S.%f')
        if not os.path.isdir("{0}/{1}".format(tmp_dir, now)):
            try:
                subprocess.check_call("mkdir {0}/{1}".format(tmp_dir, now), shell=True)
                # If another script beats us to it, this we'll fail so we'll try again then
            except:
                continue
            break

    if "out_dir_group" in settings:
        base_output_dir = "{0}/{1}/{2}".format(out_dir, settings["out_dir_group"], now)
    else:
        base_output_dir = "{0}/{1}".format(out_dir, now)
    base_output_dir = base_output_dir + "_" + config_file.split("/")[-1].split(".")[0]
    base_tmp_dir = "{0}/{1}".format(tmp_dir, now)

    # Save config file and current Git log for reproducibility.
    save_state(config_file, base_output_dir)

    if not "overlap_loci" in settings["selection_basis"]:
	dispatch_all_loci(settings, max_cores, base_output_dir, base_tmp_dir)
	
    else:
	logging.info("'overlap_loci' selection basis.")
	
	g = os.path.basename(settings["selection_basis"]["overlap_loci"]).split('.')
	g = '.'.join(g[0:len(g)-1]).replace('.','_')
	
	logging.info("Starting analysis for {}.".format(os.path.basename(settings["selection_basis"]["overlap_loci"])))
	
	def initialize_file(f, header):
	    if os.path.isfile(f):
	        raise Exception("Output file already exists: {}".format(f))
	    else: 
                with open(f, "w") as w:
                    w.write(header)
	
	# Single output file for all loci
	# Write headers
        if "finemap" in settings["methods"]:
	    out = "{0}/{1}_finemap_clpp_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tclpp\t-log_gwas_pval\t-log_eqtl_pval\tbase_gwas_file\tclpp_mod\n")
		
        if "coloc" in settings["methods"]:
	    out = "{0}/{1}_coloc_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tclpp_h0\tclpp_h1\tclpp_h2\tclpp_h3\tclpp_h4\tbase_gwas_file\n")

        if "ecaviar" in settings["methods"]:
	    out = "{0}/{1}_ecaviar_clpp_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tclpp\n")

        if "rtc" in settings["methods"]:
	    out = "{0}/{1}_rtc_score_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tgwas_trait\tfeature\trtc_score\tbase_gwas_file\n")

        if "caviarbf" in settings["methods"]:
	    out = "{0}/{1}_caviarbf_clpp_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tclpp\n")

        if "twas" in settings["methods"]:
	    out = "{0}/{1}_twas_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tfeature\tn_snps\tgwas_trait\tbase_gwas_file\ttwas_log_pval\ttwas_perm_log_pval\n")
        
        if "smr" in settings["methods"]:
	    out = "{0}/{1}_smr_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tgwas_trait\tfeature\tnum_sites\tbase_gwas_file\tsmr_neg_log_pval\theidi_pval\n")

        if "gsmr" in settings["methods"]:
	    out = "{0}/{1}_gsmr_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tgwas_trait\tfeature\tnum_sites\tbase_gwas_file\tsmr_neg_log_pval\n")

        if "metaxcan" in settings["methods"]:
	    out = "{0}/{1}_metaxcan_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tgwas_trait\tbase_gwas_file\ttwas_log_pval\n")

        if "baseline" in settings["methods"]:
	    out = "{0}/{1}_baseline_status.txt".format(base_output_dir, g)
	    initialize_file(out, "ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tbase_gwas_file\tbaseline_pval\tbaseline_pval2\tbaseline_pval3\tbaseline_pval4\tbaseline_pval5\n")

        if "ensemble" in settings["methods"]:
	    out = "{0}/{1}_ensemble_status.txt".format(base_output_dir, g)
	    intialize_file(out, "ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tbase_gwas_file\tensemble_score\n")
	
	# get total number of tests (wc of file)
        with open(settings["selection_basis"]["overlap_loci"], 'r') as f:
            for i, l in enumerate(f):
                pass
        num_tests = i + 1
	
	bar = Bar('Processing\n', max=num_tests)
	
        def update_bar(result):
	    bar.next()
	
	pool = Pool(max_cores)
	
	# iterate over loci
	with open(settings["selection_basis"]["overlap_loci"], 'r') as f:
	    header = f.readline().strip().split("\t")
	    gwas_index = header.index("gwas_file")
	    eqtl_index = header.index("eqtl_file")
	    chrom_index = header.index("chr")
	    snp_pos_index = header.index("snp_pos")
	    feature_index = header.index("feature")
	    if "trait" in header:
	        trait_index = header.index("trait")
	    else:
		trait_index = header.index("gwas_file")

	    for line in f:
		eqtl_file = line.strip().split("\t")[eqtl_index]
		gwas_file = line.strip().split("\t")[gwas_index]
		chrom = line.strip().split("\t")[chrom_index]
		snp_pos = line.strip().split("\t")[snp_pos_index]
		trait = os.path.basename(line.strip().split("\t")[trait_index])
		feature = line.strip().split("\t")[feature_index]
		
		# format SNP
		this_snp = tuple([chrom, snp_pos])
		
		# skip if it's not autosomal
		if "chr" in str(this_snp[0]):
                    try:
                        c = int(this_snp[0][3:])
                    except:
                        continue
                else:
                    try:
                        c = int(this_snp[0])
                    except:
                        continue
		
		ready_snp = SNP.SNP(this_snp) 
		pool.apply_async(analyze_snp_wrapper, 
				 args=(gwas_file, 
				       eqtl_file, 
				       ready_snp, 
				       settings, 
				       base_output_dir, 
				       base_tmp_dir, 
				       trait,
				       feature), 
				 kwds=dict(restrict_gene=-1), 
				 callback=update_bar)
	    pool.close()
	    pool.join()

	    # Clean up after ourselves
	    if not settings["debug"]:
	        subprocess.call("rm -r {0} 2> /dev/null".format(base_tmp_dir), shell=True)	
	
def dispatch_all_loci(settings, max_cores, base_output_dir, base_tmp_dir):
	
    # TODO: make this cleaner? 
	
    gwas_files = sorted([f for f in settings["gwas_experiments"]])
    eqtl_files = sorted([f for f in settings["eqtl_experiments"]])

    # For each GWAS experiment:
    for gwas_file in gwas_files:
		     
	logging.info("Initialize anlaysis for {}.".format(os.path.basename(gwas_file)))
        
        gwas_suffix = gwas_file.split("/")[-1].replace(".", "_")
        
        # Write header of output file for FINEMAP
        if "finemap" in settings["methods"]:
            with open("{0}/{1}_finemap_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tclpp\t-log_gwas_pval\t-log_eqtl_pval\tbase_gwas_file\tclpp_mod\n")

        # Write COLOC results to the desired file.
        if "coloc" in settings["methods"]:
            with open("{0}/{1}_coloc_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tclpp_h0\tclpp_h1\tclpp_h2\tclpp_h3\tclpp_h4\tbase_gwas_file\n")

        if "ecaviar" in settings["methods"]:
            with open("{0}/{1}_ecaviar_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tclpp\n")

        if "rtc" in settings["methods"]:
            with open("{0}/{1}_rtc_score_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\trtc_score\tbase_gwas_file\n")

        if "caviarbf" in settings["methods"]:
            with open("{0}/{1}_caviarbf_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tclpp\n")

        if "twas" in settings["methods"]:
            with open("{0}/{1}_twas_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tfeature\tn_snps\tgwas_trait\tbase_gwas_file\ttwas_log_pval\ttwas_perm_log_pval\n")
        
        if "smr" in settings["methods"]:
            with open("{0}/{1}_smr_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tnum_sites\tbase_gwas_file\tsmr_neg_log_pval\theidi_pval\n")

        if "gsmr" in settings["methods"]:
            with open("{0}/{1}_gsmr_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tnum_sites\tbase_gwas_file\tsmr_neg_log_pval\n")

        if "metaxcan" in settings["methods"]:
            with open("{0}/{1}_metaxcan_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tgwas_trait\tbase_gwas_file\ttwas_log_pval\n")

        if "baseline" in settings["methods"]:
            with open("{0}/{1}_baseline_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tbase_gwas_file\tbaseline_pval\tbaseline_pval2\tbaseline_pval3\tbaseline_pval4\tbaseline_pval5\n")

        if "ensemble" in settings["methods"]:
            with open("{0}/{1}_ensemble_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
                w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tbase_gwas_file\tensemble_score\n")
        
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

	    logging.info("Starting analysis for {}.".format(trait))

            gwas_snp_list = []
            # Get a list of which SNPs we should test in this GWAS.

            # Selection basis options:
            #   - gwas: SNPs significant in GWAS will be tested (at specified threshold)
            #   - eqtl: SNPs significant in eQTLs will be tested (at specified threshold)
            #   - both: SNPs significant in GWAS or eQTL will be tested
            # To run for the entire genome, specify "eQTL" and set the pvalue cutoff to 1.

            if "gwas" in settings["selection_basis"] or "both" in settings["selection_basis"]:
                gwas_snp_list.extend(preprocess.select_test_snps_by_gwas(gwas_file, settings['selection_thresholds']["gwas"], trait, settings))

            if "snps_from_list" in settings["selection_basis"]:
                if len(settings["selection_basis"]["snps_from_list"]) > 0:
                    gwas_snp_list.extend(preprocess.select_snps_from_list(settings["selection_basis"]["snps_from_list"]))
                else:
                    gwas_snp_list.extend(preprocess.select_snps_from_list(settings["gwas_experiments"][gwas_file]["snp_list_file"]))

            # For each eQTL experiment:
            for eqtl_file in eqtl_files:

                logging.info(eqtl_file)

                eqtl_snp_list = []
		if "eqtl" in settings["selection_basis"] or "both" in settings["selection_basis"]:
                    # If a "selection subset" is specified for the eQTL experiment, then genes will
                    # only be tested if they are in this subset.
                    if "selection_subset" in settings['eqtl_experiments'][eqtl_file]:
                        eqtl_subset = settings['eqtl_experiments'][eqtl_file]["selection_subset"]
                    else:
                        eqtl_subset = -1
                    eqtl_snp_list.extend(preprocess.select_test_snps_by_eqtl(eqtl_file, settings, eqtl_subset))

                snp_list = eqtl_snp_list + gwas_snp_list
                logging.info("Testing {2} SNPs ({0} GWAS hits and {1} eQTL hits).".format(len(gwas_snp_list), len(eqtl_snp_list), len(snp_list)))
		
		num_tests = len(eqtl_snp_list) + len(gwas_snp_list)

		bar = Bar('Processing', max=num_tests)

		def update_bar(result):
		    bar.next()
		
		feature = None
		
                # Run key SNPs in parallel
                pool = Pool(max_cores)
                for i in xrange(0, len(eqtl_snp_list)):
                    snp = eqtl_snp_list[i]
                    pool.apply_async(analyze_snp_wrapper, args=(gwas_file, eqtl_file, snp[0], settings, base_output_dir, base_tmp_dir, trait, feature), kwds=dict(restrict_gene=snp[1]), callback=update_bar)
                pool.close()
                pool.join()

                # Clean up after ourselves
                if not settings["debug"]:
                    subprocess.call("rm -r {0} 2> /dev/null".format(base_tmp_dir), shell=True)

                # Run GWAS SNPs separately just in case there happen to be any overlaps,
                # which could lead to a race.
                pool = Pool(max_cores)
                for i in xrange(0, len(gwas_snp_list)):
                    snp = gwas_snp_list[i]
                    pool.apply_async(analyze_snp_wrapper, args=(gwas_file, eqtl_file, snp[0], settings, base_output_dir, base_tmp_dir, trait, feature), kwds=dict(restrict_gene=snp[1]), callback=update_bar)
		pool.close()
                pool.join()

                # Clean up after ourselves
                if not settings["debug"]:
                    subprocess.call("rm -r {0} 2> /dev/null".format(base_tmp_dir), shell=True)

		bar.finish()

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
def analyze_snp_wrapper(gwas_file, eqtl_file, snp, settings, base_output_dir, base_tmp_dir, trait, feature, restrict_gene=-1):
    try:
        analyze_snp(gwas_file, eqtl_file, snp, settings, base_output_dir, base_tmp_dir, trait, feature, restrict_gene)
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        error = str(e)
        error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
        with open("{0}/ERROR_variants.txt".format(base_output_dir),"a") as a:
            a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tgeneral_error\n".format(gwas_file, eqtl_file, snp.chrom, snp.pos, restrict_gene, trait, error))
        raise Exception("Failed colocalization run.")

def analyze_snp(gwas_file, eqtl_file, snp, settings, base_output_dir, base_tmp_dir, trait, feature, restrict_gene=-1):

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

    #print genes

    # if "feature" is specified, only test this gene 
    if feature is not None:
	genes = [feature]

    # Loop through all genes now
    for gene in genes:
        # NOTE: It would be faster to just do this step once outside of this loop,
        # and then filter down to the gene of interest. Consider modifying.

        # Make sure this is a gene we actually care about
        if "selection_subset" in settings['eqtl_experiments'][eqtl_file] and gene not in settings['eqtl_experiments'][eqtl_file]["selection_subset"]:
            continue
        #print gene

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
    #subprocess.check_call('git --git-dir /users/mgloud/projects/brain_gwas/.git log >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    #subprocess.check_call('git --git-dir /users/mgloud/projects/brain_gwas/.git diff >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    #subprocess.check_call('git --git-dir /users/mgloud/projects/brain_gwas/.git branch >> {0}/git_status.txt'.format(base_output_dir), shell=True)
    #subprocess.check_call('git --git-dir /users/mgloud/projects/brain_gwas/.git status >> {0}/git_status.txt'.format(base_output_dir), shell=True)

# Use SplicePlot to generate splice plots
def splice_plot(results_file, eqtl_file, settings):
    map_file = settings["splice_plots"][eqtl_file]["map_file"]
    vcf_file = settings["splice_plots"][eqtl_file]["vcf_file"]
    sp.generate_splice_plots(results_file, eqtl_file, vcf_file, map_file, threshold = 0.01)

if __name__ == "__main__":
    main()
