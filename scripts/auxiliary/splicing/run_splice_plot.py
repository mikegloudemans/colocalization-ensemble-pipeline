#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 12/12/2017

# Given a results file of sQTL/GWAS colocalization scores, generate full splice plots
# for all loci with high probability of colocalization.

# NOTE: One major caveat of this analysis is that it will skip loci with unannotated exons;
# therefore if novel splice sites exist, we won't be able to view their splice plots.
# This occasionally happens for some of the tested loci.

# Another caveat is that the genotype used for stratifying patients is the marker variant;
# if the actual causal variant is a different one, the observed effect in the plots
# may be reduced or missing.

import sys
import subprocess

def generate_splice_plots(results_file, splice_file, vcf_file, map_file, threshold = 0.01, junction_delim=".", strict_filters=0):

    splice_suffix = splice_file.split("/")[-1].replace(".", "_")

    out_dir = "/".join(results_file.split("/")[:-1]) + "/spliceplot"
    subprocess.check_call("mkdir -p {0}".format(out_dir), shell=True)

    splice_clusters = []
    # Load data
    with open(results_file) as f:
        for line in f:
            data = line.strip().split()

            # Make sure we're only pulling information corresponding to the sQTL file of interest
            if data[1] != splice_suffix:
                continue

            score = float(data[6])
            if score < threshold:
                continue

            if float(data[7]) < strict_filters or float(data[8]) < strict_filters:
                continue

            snp = data[0]
            chrom, pos = snp.split("_")
            
            cluster = data[3].strip().split(junction_delim)[-1]
            test = (chrom, pos, cluster)
            if test not in splice_clusters:
                splice_clusters.append(test)

    for sc in splice_clusters:
        run_spliceplot(sc, splice_file, vcf_file, map_file, out_dir, junction_delim)

# Make separate function to call spliceplot for a single locus; that way someone
# can call spliceplot with a particular cluster in mind, or submit a batch that
# does a bunch in a row.

def run_spliceplot(splice_site, splice_file, vcf_file, map_file, out_dir, junction_delim="."):

    # Read header information to understand column organization of splice file
    header = subprocess.check_output("zcat {0} | head -n 1".format(splice_file, shell=True), shell=True).strip().split()
    chr_index = header.index("chr")
    pos_index = header.index("snp_pos")
    feature_index = header.index("feature")

    sc = splice_site
    print sc
    # Get all splice junctions in this cluster
    junctions = []
    first = {}
    second = {}
    text = subprocess.check_output("zcat {0} | grep {1} | grep {2}".format(splice_file, sc[2], sc[1]), shell=True)

    for line in text.strip().split("\n"):

        data = line.replace(":", ".").strip().split()

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

        # Hive plot would fail if we have more than 2 possible splicing events
        # This might be fixable by modifying the config file, but it's not worth it right now
        no_hive = ""
        if len(tj) > 2:
            no_hive = "_nohive"

        # Initialize data for splice plot in pickle format
        print "python /users/mgloud/projects/github_forks/SplicePlot/initialize_data.py chr{0}:{1} {2} --vcf {3} --gtf /users/mgloud/projects/brain_gwas/data/spliceplot/hg19.exons.gtf.gz --mf {4} --output /users/mgloud/projects/github_forks/SplicePlot/pickle_files/preplot.p".format(sc[0], sc[1], junction_text, vcf_formatted, map_file)
        subprocess.check_call("python /users/mgloud/projects/github_forks/SplicePlot/initialize_data.py chr{0}:{1} {2} --vcf {3} --gtf /users/mgloud/projects/brain_gwas/data/spliceplot/hg19.exons.gtf.gz --mf {4} --output /users/mgloud/projects/github_forks/SplicePlot/pickle_files/preplot.p".format(sc[0], sc[1], junction_text, vcf_formatted, map_file), shell=True)

        # Create plots
        print "python /users/mgloud/projects/github_forks/SplicePlot/plot.py /users/mgloud/projects/github_forks/SplicePlot/pickle_files/preplot.p pickle /users/mgloud/projects/github_forks/SplicePlot/settings_file{5} --output {0}/{1}-{2}-{3}_{4}".format(out_dir, sc[0], sc[1], sc[2], junction_text.replace(",", "_").replace(":", "-"), no_hive)
        #status = subprocess.check_output("python /users/mgloud/projects/github_forks/SplicePlot/plot.py /users/mgloud/projects/github_forks/SplicePlot/pickle_files/preplot.p pickle /users/mgloud/projects/github_forks/SplicePlot/settings_file{5} --output {0}/{1}-{2}-{3}_{4}".format(out_dir, sc[0], sc[1], sc[2], junction_text.replace(",", "_").replace(":", "-"), no_hive), shell=True)
        status = subprocess.check_output("python /users/mgloud/projects/github_forks/SplicePlot/plot.py /users/mgloud/projects/github_forks/SplicePlot/pickle_files/preplot.p pickle /users/mgloud/projects/github_forks/SplicePlot/settings_file_scaled --output {0}/{1}-{2}-{3}_{4}".format(out_dir, sc[0], sc[1], sc[2], junction_text.replace(",", "_").replace(":", "-"), no_hive), shell=True)

        print status

        # TODO: Warning: In current form, the SplicePlot call can crash but the program will
        # still keep running. Fix this to make sure it worked, maybe by scanning the printed output for the words "Done!" or "Failed".
        # (This is hopefully fixed now)
        assert "Done!" in status



if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "Usage:\nrun_splice_plot.py [sqtl_chrom] [sqtl_position] [sqtl_target_cluster] [splice_file] [vcf_file] [map_file] [splice_site_delimiter] [optional:output_directory]"
        sys.exit()

    sc = (sys.argv[1], sys.argv[2], sys.argv[3])
    splice_file = sys.argv[4]
    vcf_file = sys.argv[5]
    map_file = sys.argv[6]

    junction_delim = sys.argv[7]

    try:
        out_dir = sys.argv[8]
    except:
        out_dir = ""

    # We can also call spliceplot from the command line, for a single site.
    run_spliceplot(sc, splice_file, vcf_file, map_file, out_dir, junction_delim)

