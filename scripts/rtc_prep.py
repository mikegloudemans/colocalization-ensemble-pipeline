# Functions required to prepare for running RTC for the first time on a dataset
# Author: Mike Gloudemans
# Date created: 6/2/2020

import glob
import pandas as pd
import pyarrow
import hashlib
import os
import subprocess

# Files with the info we need
gtex_v8_vcf_file = "data/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"
gtex_v8_eqtl_phenos_file = "data/gtex-reformatted/eqtl.tpm.txt"
gtex_v8_sqtl_phenos_file = "data/gtex/GTEx_Analysis_2017-06-05_v8/rna_seq/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.ft"
gtex_v8_tissue_mapping_file = "data/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"

# Files to create
rtc_index_dir = "data/rtc-gene-index"

def index_gtex_v8():

    # May first require using R to transfer from feather to a regular table, at least temporarily

    eqtl_tissues = glob.glob("data/eqtls/gtex_v8/*")
    for tissue in eqtl_tissues:
        
        short_tissue = tissue.split(".")[0].split("/")[-1]
        print short_tissue

        subset_ids = get_samp_ids(short_tissue, gtex_v8_tissue_mapping_file)

        index_pheno_file(gtex_v8_eqtl_phenos_file, subset_ids, short_tissue, rtc_index_dir + "/eqtls")
        index_geno_file(gtex_v8_vcf_file, subset_ids, short_tissue, rtc_index_dir + "/genos")

# Return a list of the GTEx samples that belong to this particular tissue
def get_samp_ids(short_tissue, gtex_v8_tissue_mapping_file):
    # Load in map of IDs
    id_map = pd.read_csv(gtex_v8_tissue_mapping_file, sep="\t")
    id_map['SMTSD'] = id_map['SMTSD'].astype(str).apply(lambda x: x.replace(" - ", "_")).apply(lambda x: x.replace(" ", "_")).apply(lambda x: x.replace("(", "")).apply(lambda x: x.replace(")", ""))

    id_map = id_map[(id_map['SMTSD'] == short_tissue).values]
    tissue_samples = list(id_map['SAMPID'])
    
    return tissue_samples


def index_pheno_file(source_file, individual_subset, tissue, out_folder):

    data = pd.read_csv(source_file, sep="\t")
 
    # Get rid of any individuals that for some reason don't have 
    # phenotypes listed (maybe there was no RNA-seq)
    sub_sub = ["gene_id"] + [inds for inds in individual_subset if inds in data.columns.values]

    tissue_data = data[sub_sub]

    if not os.path.exists("{0}".format(out_folder)):
        os.mkdir("{0}".format(out_folder))

    with open("{0}/header.txt".format(out_folder), "w") as w:
        w.write("\t".join(tissue_data.columns.values))

    # For each line in the file
    for i in range(tissue_data.shape[0]):

        gene = data['gene_id'][i]
        hsh = hashlib.sha1(gene).hexdigest()
        if not os.path.exists("{0}/{1}".format(out_folder, hsh[:2])):
            os.mkdir("{0}/{1}".format(out_folder, hsh[:2]))
        if not os.path.exists("{0}/{1}/{2}".format(out_folder, hsh[:2], hsh[2:4])):
            os.mkdir("{0}/{1}/{2}".format(out_folder, hsh[:2], hsh[2:4]))
        with open("{0}/{1}/{2}/{3}_{4}.txt".format(out_folder, hsh[:2], hsh[2:4], gene, tissue), "w") as w:
            w.write("\t".join([str(s) for s in list(tissue_data.iloc[i,:])]))


def index_geno_file(source_file, individual_subset, tissue, out_folder):

    # Load the original genotypes file, contained in source_file
    data = pd.read_csv(gtex_v8_vcf_file, sep="\t", skiprows=3385)

    # Subset it to the individuals specified in individual_subset
    subjects = ["#CHROM", "POS", "ID", "REF", "ALT"] + ["-".join(ind.split("-")[:2]) for ind in individual_subset if "-".join(ind.split("-")[:2]) in data.columns.values]
    tissue_data = data[subjects]

    if not os.path.exists("{0}".format(out_folder)):
        os.mkdir("{0}".format(out_folder))

    # Write it to out_folder/out_file for subsequent access
    with open("{0}/{1}.txt".format(out_folder, tissue), "w") as w:
        tissue_data.to_csv(w, sep="\t", index=False)

    subprocess.check_call("bgzip -f {0}/{1}.txt".format(out_folder, tissue), shell=True)
    subprocess.check_call("tabix -S 1 -s 1 -b 2 -e 2 {0}/{1}.txt.gz".format(out_folder, tissue), shell=True)

def main():
    index_gtex_v8()

if __name__ == "__main__":
    main()
