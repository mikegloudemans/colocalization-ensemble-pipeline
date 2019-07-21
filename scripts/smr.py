# Author: Mike Gloudemans
# Author: Abhiram Rao

import subprocess
from scipy import stats
from shutil import copyfile
import sys
if sys.version_info[0] < 3: 
   from StringIO import StringIO
else:
   from io import StringIO
import pandas as pd
import math
import os
import numpy as np

# Integration of code written by Ram
def run_smr(locus, window=1000000):
    
    subprocess.call("mkdir -p {4}/smr/{0}/{1}_{2}/{3}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.tmpdir), shell=True)

    dat = locus.data.copy()

    # We need beta and se for GWAS, or else this analysis won't work
    if "beta_gwas" not in dat or "se_gwas" not in dat:
        raise Exception("Need beta and SE for GWAS in order to run SMR.")

    # Get rsids
    # TODO later: Move this to an earlier point so it doesn't have to be repeated
    # for plotting or anything
    
    # If both GWAS and eQTL have rsids, they'll have been renamed
    if "rsid_eqtl" in list(dat.columns.values):
        dat['rsid'] = dat['rsid_eqtl']

    if "rsid" not in list(dat.columns.values):
        if "rsid_index_file" in locus.settings:

            # First, extract nearby variants using tabix
            stream = StringIO(subprocess.check_output("tabix {0} {1}:{2}-{3}".format(locus.settings["rsid_index_file"], locus.data['chr_gwas'][0], np.min(np.array(locus.data["snp_pos"])), np.max(np.array(locus.data["snp_pos"]))), shell=True))

            # For readability, load the header too
            # Load with pandas
            dbsnp = pd.read_csv(stream, sep="\t", header=None).iloc[:,:3]
            dbsnp = dbsnp.rename({0: "chr_gwas", 1: "snp_pos", 2:"rsid"}, axis="columns")
           
            dat = pd.merge(dbsnp, dat, left_on=["chr_gwas", "snp_pos"], right_on=["chr_gwas", "snp_pos"])

        else:
            raise Exception("No rsids specified; skipping locus.")

    # See if we need to infer effect_af from the reference VCF
    if "effect_af_gwas" not in dat.columns.values or "effect_af_eqtl" not in dat.columns.values:
        # Code for this is currently in the COLOC script
        import coloc
        dat = coloc.get_mafs(locus, dat, window)
        if "effect_af_eqtl" not in dat:
            dat["effect_af_eqtl"] = dat["ref_af"]
        if "effect_af_gwas" not in dat:
            dat["effect_af_gwas"] = dat["ref_af"]

    smr_location = '/users/raoa/coloc_comparison/SMR_v1.01/smr_Linux'
    ld_reference = '/users/raoa/coloc_comparison/ld_reference'
    eqtl_pval_threshold = 1 # p-value threshold so that we run SMR at every test
    maf_threshold = 0.01
    cis_wind = 10000 # large cis window to include all SNPs
    diff_freq = 1 # allowable allele freq difference between datasets
    diff_freq_prop = 1 # all SNPs are allowed to have allele freq differences of at least diff.freq

    # probably don't actually need to do these checks
    '''if min(dat.chr) != max(dat.chr):
	print("Multiple chromosomes in eQTL input file")

    if min(dat.chr) != max(dat.chr):
	print("Multiple chromosomes in eQTL input file")

    if len(dat.feature.unique()) > 1:
	print("There are multiple genes in the eqtl input file")'''

    # Run SMR and get pvalue

    # inferred info
    gene = dat.gene.unique()[0]
    chrom = min(dat.chr_gwas)
    tss = np.round(np.mean(dat.snp_pos)) # Probably not exactly correct, but good enough for now
    
    ## create esd
    esd = pd.DataFrame({'Chr': dat.chr_eqtl, 'SNP': dat.rsid, 'Bp': dat.snp_pos, 'A1': dat.alt_eqtl, 'A2': dat.ref_eqtl, 'Freq': dat.effect_af_eqtl, 'Beta': dat.beta_eqtl, 'se': dat.se_eqtl, 'p': dat.pvalue_eqtl})
    esd = esd[['Chr', 'SNP', 'Bp', 'A1', 'A2', 'Freq', 'Beta', 'se', 'p']]
    esd.to_csv('{6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.esd'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), sep="\t", index = False)

    ## create flist file
    flist = pd.DataFrame({'Chr': chrom, 'ProbeID': gene, 'GeneticDistance': '0', 'ProbeBp': tss, 'Gene': gene, 'Orientation': '+', 'PathOfEsd': '{6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.esd'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)}, index=[0])
    flist = flist[['Chr', 'ProbeID', 'GeneticDistance', 'ProbeBp', 'Gene', 'Orientation', 'PathOfEsd', ]]
    flist.to_csv('{6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.flist'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), sep="\t", index = False)

    ## run SMR to create binary files
    os.system('{7} --eqtl-flist {6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.flist --make-besd --out {6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.besd > /dev/null'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir, smr_location))

    ## run SMR with corresponding gwas file 

    # TODO: TODO: Warning! Currently we're just guessing the value of N! It would be
    # wise to get the actual values of N for these studies so that we 
    # can use them for COLOC, SMR, any other tools where N is needed.
    if "n_cases" in dat and "n_controls" in dat:
        n = dat.n_cases + dat.n_controls
    else:
        n = 10000
    dat['n'] = n

    gwas = dat[['rsid','alt_gwas','ref_gwas','effect_af_gwas','beta_gwas','se_gwas','pvalue_gwas','n']]
    gwas.columns = ['SNP','A1','A2','freq','b','se','p','n']
    gwas.to_csv("{6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.gwas".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir, smr_location), sep="\t", index = False)
    os.system('/users/raoa/coloc_comparison/SMR_v1.01/smr_Linux --bfile /users/raoa/coloc_comparison/ld_reference/EUR1KG_chr{1}_hg19_rsid --gwas-summary {6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.gwas --beqtl-summary {6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.besd --maf {7} --thread-num 1 --diff-freq {8} --diff-freq-prop {9} --peqtl-smr {10} --cis-wind {11} --out {6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.out > /dev/null'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir, maf_threshold, diff_freq, diff_freq_prop, eqtl_pval_threshold, cis_wind))

    # Now parse results
    with open('{6}/smr/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_smr.out.smr'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)) as f:
        f.readline()
        data = f.readline().strip().split()
        pval = float(data[18])

    # Write results to the output file
    with open("{0}/{1}_smr_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, locus.data.shape[0], locus.gwas_suffix, -1*math.log10(pval), data[19]))

    return pval

