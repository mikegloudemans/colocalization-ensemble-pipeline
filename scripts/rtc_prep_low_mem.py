# Functions required to prepare for running RTC for the first time on a dataset
# Author: Mike Gloudemans
# Date created: 6/2/2020

import glob
import pandas as pd
import hashlib
import os
import subprocess
import gzip

# Files with the info we need
gtex_v8_vcf_file = "data/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"
gtex_v8_eqtl_phenos_file = "data/gtex-reformatted/eqtl.tpm.txt"
gtex_v8_sqtl_phenos_file = "data/gtex/GTEx_Analysis_2017-06-05_v8/rna_seq/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.ft"
gtex_v8_tissue_mapping_file = "data/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"

# Files to create
rtc_index_dir = "data/rtc-gene-index"

def index_gtex_v8():

    # May first require using R to transfer from feather to a regular table, at least temporarily

    eqtl_tissues = [g.split(".")[0].split("/")[-1] for g in glob.glob("data/eqtls/gtex_v8/*")]

    index_pheno_file(gtex_v8_eqtl_phenos_file, eqtl_tissues, rtc_index_dir + "/expression")
    index_geno_file(gtex_v8_vcf_file, eqtl_tissues, rtc_index_dir + "/genotypes")

# Return a list of the GTEx samples that belong to this particular tissue
def get_samp_ids(short_tissue, gtex_v8_tissue_mapping_file):
    # Load in map of IDs
    id_map = pd.read_csv(gtex_v8_tissue_mapping_file, sep="\t")
    id_map['SMTSD'] = id_map['SMTSD'].astype(str).apply(lambda x: x.replace(" - ", "_")).apply(lambda x: x.replace(" ", "_")).apply(lambda x: x.replace("(", "")).apply(lambda x: x.replace(")", ""))

    id_map = id_map[(id_map['SMTSD'] == short_tissue).values]
    tissue_samples = list(id_map['SAMPID'])
    
    return tissue_samples


def get_samp_id_dict(all_tissues, gtex_v8_tissue_mapping_file):
    samp_id_dict = {}
    for tissue in all_tissues:
        samp_id_dict[tissue] = get_samp_ids(tissue, gtex_v8_tissue_mapping_file)
    return samp_id_dict


def index_pheno_file(source_file, tissue_list, out_folder):
    if not os.path.exists("{0}".format(out_folder)):
        os.mkdir("{0}".format(out_folder))

    # Figure out which columns need to be saved for each tissue
    tissue_ids = get_samp_id_dict(tissue_list, gtex_v8_tissue_mapping_file)
    tissue_columns = {}
    with open(source_file) as f:
        header = f.readline().strip().split()

#chr    start   end     gene    length  strand

        for tissue in tissue_ids:
            tissue_columns[tissue] = [i for i,e in enumerate(header) if e in ["gene_id"] + tissue_ids[tissue]]

            # Write the header file while we're at it
            if not os.path.exists("{0}/{1}".format(out_folder, tissue)):
                os.mkdir("{0}/{1}".format(out_folder, tissue))
            with open("{0}/{1}/header.txt".format(out_folder, tissue), "w") as w:
                w.write("\t".join([header[i] for i in tissue_columns[tissue]]))

        if not os.path.exists("{0}".format(out_folder)):
            os.mkdir("{0}".format(out_folder))

        # Now go through the whole file
        for line in f:
            data = line.strip().split()
            gene = data[0]
            hsh = hashlib.sha1(gene).hexdigest()
            print gene
            # Tackle one tissue at a time
            for tissue in tissue_columns:
                if not os.path.exists("{0}/{1}/{2}".format(out_folder, tissue, hsh[:2])):
                    os.mkdir("{0}/{1}/{2}".format(out_folder, tissue, hsh[:2]))
                with open("{0}/{1}/{2}/{3}.txt".format(out_folder, tissue, hsh[:2], gene), "w") as w:
                    tissue_data = [data[i] for i in tissue_columns[tissue]]
                    w.write("\t".join([str(s) for s in tissue_data]))

def index_geno_file(source_file, tissue_list, out_folder):

    if not os.path.exists("{0}".format(out_folder)):
        os.mkdir("{0}".format(out_folder))

    # Figure out which columns need to be saved for each tissue
    tissue_ids = get_samp_id_dict(tissue_list, gtex_v8_tissue_mapping_file)
    for tissue in tissue_ids:
        tissue_ids[tissue] = ["-".join(ti.split("-")[:2]) for ti in tissue_ids[tissue]]
    
    tissue_columns = {}
    with gzip.open(source_file) as f:
        header = f.readline().strip()
        while header.startswith("##"):
            header = f.readline().strip()

        header = header.split()
        # Write the header first
        for tissue in tissue_ids:
            tissue_columns[tissue] = [i for i,e in enumerate(header) if e in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + tissue_ids[tissue]]

            # Write the header file while we're at it
            subprocess.check_call("rm -f {0}/{1}.txt".format(out_folder, tissue), shell=True)
            with open("{0}/{1}.txt".format(out_folder, tissue), "w") as w:
                w.write("##fileformat=VCFv4.1\n")
                w.write("\t".join([header[i] for i in tissue_columns[tissue]]) + "\n")

        # Now go through the whole file
        for line in f:
            data = line.strip().split()
            print data[:4]

            for tissue in tissue_columns:
                # Write it to out_folder/out_file for subsequent access
                with open("{0}/{1}.txt".format(out_folder, tissue), "a") as a:
                    a.write("\t".join([data[s] for s in tissue_columns[tissue]]) + "\n")

    for tissue in tissue_columns:
        subprocess.check_call("bgzip -f {0}/{1}.txt".format(out_folder, tissue), shell=True)
        subprocess.check_call("tabix -S 1 -s 1 -b 2 -e 2 {0}/{1}.txt.gz".format(out_folder, tissue), shell=True)

def main():
    index_gtex_v8()

if __name__ == "__main__":
    main()
