#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 12/12/2017

# Given a results file of sQTL/GWAS colocalization scores, generate full splice plots
# for all loci with high probability of colocalization.

# NOTE: One major caveat of this analysis is that it will skip loci with unannotated exons;
# therefore if novel splice sites exist, we won't be able to view their splice plots.

import sys
import subprocess

def generate_splice_plots(results_file, splice_file, vcf_file, map_file, threshold = 0.01, junction_delim=":"):

    out_dir = "/".join(results_file.split("/")[:-1]) + "/spliceplot"
    subprocess.check_call("mkdir -p {0}".format(out_dir), shell=True)

    splice_clusters = []
    # Load data
    with open(results_file) as f:
        for line in f:
            data = line.strip().split()
            score = float(data[6])
            if score < threshold:
                continue

            snp = data[0]
            chrom, pos = snp.split("_")
            
            cluster = data[3].strip().split(junction_delim)[-1]
            test = (chrom, pos, cluster)
            if test not in splice_clusters:
                splice_clusters.append(test)

    # Read header information to understand column organization of splice file
    header = subprocess.check_output("zcat {0}".format(splice_file, shell=True), shell=True).strip().split()
    chr_index = header.index("chr")
    pos_index = header.index("snp_pos")
    feature_index = header.index("feature")
    # For each cluster passing the threshold
    for sc in splice_clusters:
        # Get all splice junctions in this cluster
        junctions = []
        first = {}
        second = {}
        text = subprocess.check_output("zcat {0} | grep {1} | grep {2}".format(splice_file, sc[2], sc[1]), shell=True)
        for line in text.strip().split("\n"):
            print line
            data = line.strip().split()

            # Just in case grep matched a different column than the one
            # we wanted
            if data[chr_index] != sc[0] or data[pos_index] != sc[1]:
                continue

            # Get data on splice junction
            junc_data = data[feature_index].strip().split(junction_delim)
            junc_start = junc_data[1]
            junc_end = junc_data[2]

            # Check to make sure the splice junction actually corresponds with an exon in our annotation
            valid = int(subprocess.check_output("zcat /users/mgloud/projects/brain_gwas/data/spliceplot/hg19.exons.gtf.gz | awk '{{if ($5=={0}) print $0}}' | wc -l".format(junc_start),shell=True).strip()) > 0 and int(subprocess.check_output("zcat /users/mgloud/projects/brain_gwas/data/spliceplot/hg19.exons.gtf.gz | awk '{{if ($4=={0}) print $0}}' | wc -l".format(junc_end),shell=True).strip()) > 0
            
            if not valid:
                continue

            first[junc_start] = first.get(junc_start, 0) + 1
            second[junc_end] = second.get(junc_end, 0) + 1
            junctions.append((junc_start, junc_end))

        test_junctions = []
        # Find all sets of splice junctions sharing a common
        # donor or common acceptor site position
        # Shared first position
        for f in first.keys():
            if first[f] > 1:
                test_junctions.append([j for j in junctions if j[0] == f])

        # Shared second position
        for s in second.keys():
            if second[s] > 1:
                test_junctions.append([j for j in junctions if j[1] == s])

        # Run SplicePlot for each of these sets, outputting
        # results to a splice plot directory in the output folder
        print sc, test_junctions
        vcf_formatted = vcf_file.format(sc[0])

        for tj in test_junctions:

            junction_text = []
            for j in tj:
                junction_text.append("chr{0}:{1}-{2}".format(sc[0], j[0], j[1]))

            junction_text = ",".join(junction_text)

            print "python /users/mgloud/software/SplicePlot/initialize_data.py chr{0}:{1} {2} --vcf {3} --gtf /users/mgloud/projects/brain_gwas/data/spliceplot/hg19.exons.gtf.gz --mf {3} --output /users/mgloud/software/SplicePlot/pickle_files/preplot.p".format(sc[0], sc[1], junction_text, vcf_formatted, map_file)

            # Initialize data for splice plot in pickle format
            subprocess.check_call("python /users/mgloud/software/SplicePlot/initialize_data.py chr{0}:{1} {2} --vcf {3} --gtf /users/mgloud/projects/brain_gwas/data/spliceplot/hg19.exons.gtf.gz --mf {4} --output /users/mgloud/software/SplicePlot/pickle_files/preplot.p".format(sc[0], sc[1], junction_text, vcf_formatted, map_file), shell=True)

            print "python /users/mgloud/software/SplicePlot/plot.py /users/mgloud/software/SplicePlot/pickle_files/preplot.p pickle /users/mgloud/software/SplicePlot/settings_file --output {0}/{1}-{2}-{3}__{4}".format(out_dir, sc[0], sc[1], sc[2], junction_text.replace(",", "_").replace(":", "-"))

            # Create plots
            status = subprocess.check_call("python /users/mgloud/software/SplicePlot/plot.py /users/mgloud/software/SplicePlot/pickle_files/preplot.p pickle /users/mgloud/software/SplicePlot/settings_file --output {0}/{1}-{2}-{3}__{4}".format(out_dir, sc[0], sc[1], sc[2], junction_text.replace(",", "_").replace(":", "-")), shell=True)

            # TODO: Warning: In current form, the SplicePlot call can crash but the program will
            # still keep running. Fix this to make sure it worked, maybe by scanning the printed output for the words "Done!" or "Failed"

            print status
