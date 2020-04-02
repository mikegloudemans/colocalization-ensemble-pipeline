    # For each GWAS experiment:
    for gwas_file in gwas_files:
        
        #gwas_suffix = gwas_file.split("/")[-1].replace(".", "_")
        
        # # Write header of output file for FINEMAP
        # if "finemap" in settings["methods"]:
        #     with open("{0}/{1}_finemap_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tclpp\t-log_gwas_pval\t-log_eqtl_pval\tbase_gwas_file\tclpp_mod\n")

        # # Write COLOC results to the desired file.
        # if "coloc" in settings["methods"]:
        #     with open("{0}/{1}_coloc_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tclpp_h0\tclpp_h1\tclpp_h2\tclpp_h3\tclpp_h4\tbase_gwas_file\n")

        # if "ecaviar" in settings["methods"]:
        #     with open("{0}/{1}_ecaviar_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tclpp\n")

        # if "rtc" in settings["methods"]:
        #     with open("{0}/{1}_rtc_score_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\trtc_score\tbase_gwas_file\n")

        # if "caviarbf" in settings["methods"]:
        #     with open("{0}/{1}_caviarbf_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tclpp\n")

        # if "twas" in settings["methods"]:
        #     with open("{0}/{1}_twas_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tfeature\tn_snps\tgwas_trait\tbase_gwas_file\ttwas_log_pval\ttwas_perm_log_pval\n")
        
        # if "smr" in settings["methods"]:
        #     with open("{0}/{1}_smr_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tnum_sites\tbase_gwas_file\tsmr_neg_log_pval\theidi_pval\n")

        # if "gsmr" in settings["methods"]:
        #     with open("{0}/{1}_gsmr_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tnum_sites\tbase_gwas_file\tsmr_neg_log_pval\n")

        # if "metaxcan" in settings["methods"]:
        #     with open("{0}/{1}_metaxcan_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tfeature\tconditional_level\tnum_sites\tgwas_trait\tbase_gwas_file\ttwas_log_pval\n")

        # if "baseline" in settings["methods"]:
        #     with open("{0}/{1}_baseline_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tbase_gwas_file\tbaseline_pval\tbaseline_pval2\tbaseline_pval3\tbaseline_pval4\tbaseline_pval5\n")

        # if "ensemble" in settings["methods"]:
        #     with open("{0}/{1}_ensemble_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        #         w.write("ref_snp\teqtl_file\tgwas_trait\tfeature\tn_snps\tbase_gwas_file\tensemble_score\n")
        
        # # Get list of traits measured in this GWAS
        # traits = set([])

        # # Are there multiple traits in this GWAS?
        # if "traits" not in settings["gwas_experiments"][gwas_file]:
        #     with gzip.open(gwas_file) as f:
        #         header = f.readline().strip().split("\t")
        #         if "trait" in header:
        #             trait_index = header.index("trait")
        #             for line in f:
        #                 traits.add(line.strip().split("\t")[trait_index])
        #         else:
        #             traits.add(gwas_file.split("/")[-1])
        #     traits = list(traits)

        # # Subset down to traits of interest, if specified
        # if "traits" in settings["gwas_experiments"][gwas_file]:
        #     traits = settings["gwas_experiments"][gwas_file]["traits"]

        # assert len(traits) != 0

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
                gwas_snp_list.extend(preprocess.select_test_snps_by_gwas(gwas_file, settings['selection_thresholds']["gwas"], trait, settings))

            if settings["selection_basis"] == "snps_from_list":
                if "snp_list_file" in settings:
                    gwas_snp_list.extend(preprocess.select_snps_from_list(settings["snp_list_file"]))
                else:
                    gwas_snp_list.extend(preprocess.select_snps_from_list(settings["gwas_experiments"][gwas_file]["snp_list_file"]))

            # For each eQTL experiment:
            for eqtl_file in eqtl_files:

                print eqtl_file

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
		
		num_tests = len(eqtl_snp_list) + len(gwas_snp_list)

		bar = Bar('Processing\n', max=num_tests)

		def update_bar(result):
		    bar.next()
			
                # Run key SNPs in parallel
                pool = Pool(max_cores)
                for i in xrange(0, len(eqtl_snp_list)):
                    snp = eqtl_snp_list[i]
                    pool.apply_async(analyze_snp_wrapper, args=(gwas_file, eqtl_file, snp[0], settings, base_output_dir, base_tmp_dir, trait), kwds=dict(restrict_gene=snp[1]), callback=update_bar)
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
                    pool.apply_async(analyze_snp_wrapper, args=(gwas_file, eqtl_file, snp[0], settings, base_output_dir, base_tmp_dir, trait), kwds=dict(restrict_gene=snp[1]), callback=update_bar)
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

