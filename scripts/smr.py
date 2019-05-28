# Author: Mike Gloudemans

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

def run_smr(locus):
    
    combined = locus.data.copy()

    # Run SMR and get pvalue
    pval = 1

    # Write results to the output file

    with open("{0}/{1}_smr_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
        pass
        #a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.trait, locus.gene, locus.data.shape[0], locus.gwas_suffix, -1*math.log10(pval), -1*math.log10(pval2), -1*math.log10(pval3+1e-150), -1*math.log10(pval4), -1*math.log10(pval5)))

    return pval

#/usr/bin/env python
'''
eqtl_file = 'eqtl_sumstats{0}.txt.gz'.format(sim_num)
gwas_file = 'gwas_sumstats{0}.txt.gz'.format(sim_num)

# options required - config?
smr_location = '/users/raoa/coloc_comparison/SMR_v1.01/smr_Linux'
eqtl_data = '/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/eqtl'
gwas_data = '/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/gwas'
ld_reference = '/users/raoa/coloc_comparison/ld_reference'
smr_results = '/users/raoa/coloc_comparison/smr_results'
eqtl_pval_threshold = 1 # p-value threshold so that we run SMR at every test
thread_num = 10
maf_threshold = 0.01
cis_wind = 3000 # large cis window to include all SNPs
diff_freq = 1 # allowable allele freq difference between datasets
diff_freq_prop = 1 # all SNPs are allowed to have allele freq differences of at least diff.freq

# read in eqtl data
os.system('scp {0}/{1} .'.format(eqtl_data, eqtl_file))
os.system('gunzip {0}'.format(eqtl_file))
eqtl_file = eqtl_file.rstrip('.gz')
dat = pd.read_csv(eqtl_file, delimiter="\t")

# random checks
if min(dat.chr) != max(dat.chr):
	print("Multiple chromosomes in eQTL input file")

if min(dat.chr) != max(dat.chr):
	print("Multiple chromosomes in eQTL input file")

if len(dat.feature.unique()) > 1:
	print("There are multiple genes in the eqtl input file")

# inferred info
gene = dat.feature.unique()[0]
chrom = min(dat.chr)
tss = np.round(np.mean(dat.snp_pos))

## create esd
esd = pd.DataFrame({'Chr': dat.chr, 'SNP': dat.rsid, 'Bp': dat.snp_pos, 'A1': dat.alt, 'A2': dat.ref, 'Freq': dat.effect_af, 'Beta': dat.beta, 'se': dat.se, 'p': dat.pvalue})
esd = esd[['Chr', 'SNP', 'Bp', 'A1', 'A2', 'Freq', 'Beta', 'se', 'p']]
esd.to_csv('{0}{1}.esd'.format(gene, sim_num), sep="\t", index = False)

## create flist file
flist = pd.DataFrame({'Chr': chrom, 'ProbeID': gene, 'GeneticDistance': '0', 'ProbeBp': tss, 'Gene': gene, 'Orientation': '+', 'PathOfEsd': './{0}{1}.esd'.format(gene, sim_num)}, index=[0])
flist = flist[['Chr', 'ProbeID', 'GeneticDistance', 'ProbeBp', 'Gene', 'Orientation', 'PathOfEsd', ]]
flist.to_csv('{0}{1}.flist'.format(gene, sim_num), sep="\t", index = False)

## run SMR to create binary files
os.system('{0} --eqtl-flist {1}{2}.flist --make-besd --out {1}{2}'.format(smr_location, gene, sim_num))

## run SMR with corresponding gwas file
os.system('scp {0}/{1} .'.format(gwas_data, gwas_file))
os.system('gunzip {0}'.format(gwas_file))
gwas_file = gwas_file.rstrip('.gz')
gwas = pd.read_csv(gwas_file, delimiter='\t')
n = gwas.n_cases + gwas.n_controls
gwas['n'] = n
gwas = gwas[['rsid','alt','ref','effect_af','beta','se','pvalue','n']]
gwas.columns = ['SNP','A1','A2','freq','b','se','p','n']
gwas.to_csv(gwas_file, sep="\t", index = False)
os.system('{0} --bfile {1}/EUR1KG_chr{2}_hg19_rsid --gwas-summary {3} --beqtl-summary {4}{5} --maf {6} --thread-num {7} --diff-freq {8} --diff-freq-prop {9} --peqtl-smr {10} --cis-wind {11} --out {12}/sim_{5}'.format(smr_location, ld_reference, chrom, gwas_file, gene, sim_num, maf_threshold, thread_num, diff_freq, diff_freq_prop, eqtl_pval_threshold, cis_wind, smr_results))

## delete smr binaries and esd for current gene, gwas and eqtl raw files
os.system('rm {0}{1}*'.format(gene, sim_num))
os.system('rm {0}'.format(eqtl_file))
os.system('rm {0}'.format(gwas_file))
os.system('rm {0}/sim_{1}.snp_failed_freq_ck.list'.format(smr_results, sim_num))
'''
