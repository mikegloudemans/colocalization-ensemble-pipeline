#/usr/bin/env python

## Github links
# https://github.com/hakyimlab/PredictDBPipeline/wiki/Detailed_Description - this is a general description but matched to the older predictdb pipeline
# https://groups.google.com/forum/#!msg/predixcanmetaxcan/TkBxYkUpNGw/Q_mMApRtCQAJ - here's a general description of new intermediate file formats matching the new pipeline
# https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7
# https://github.com/hakyimlab/MetaXcan/wiki/Command-Line-Reference
# https://github.com/hakyimlab/MetaXcan
# uses a modified create_model.R from predictDB pipeline, also modified snp_annot_to_RDS.R and GTEx_Tissue_Wide_CV_elasticNet.R

## Steps
# clone https://github.com/hakyimlab/MetaXcan into working directory
# set working directory below, all other paths should be relative
# PredictDB_pipeline requires some modifications to run here - the modified scripts are provided with this code

import os
import subprocess
import pandas as pd
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('agg')

# arguments
sim_num = sys.argv[1]
# sim_num = 320
eqtl_file = 'eqtl_sumstats{0}.txt.gz'.format(sim_num) ## take in argument of eQTL db and covariance file?
gwas_file = 'gwas_sumstats{0}.txt.gz'.format(sim_num)
eqtl_genotypes = 'eqtl_genotypes{0}.vcf.gz'.format(sim_num)
eqtl_phenotypes = 'eqtl_phenotypes{0}.bed.gz'.format(sim_num)

# options required - config?
working_directory = '/users/raoa/coloc_comparison'
predixcan_location =  working_directory + '/MetaXcan/software'
predictdb_pipeline_location = working_directory + '/PredictDB_Pipeline_GTEx_v7' # location of PredictDB github repository
# predictdb_pipeline_scripts = predictdb_pipeline_location + '/scripts'
eqtl_data = '/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/eqtl'
gwas_data = '/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/gwas'
predixcan_results = working_directory + '/predixcan_results'
vcftools_location = working_directory + '/vcftools_0.1.13/cpp'
create_predictdb = True # required for simulations, could be substituted for downloaded dbs for real data
os.chdir(working_directory)

####### GWAS
gwas = pd.read_csv('{0}/{1}'.format(gwas_data, gwas_file), delimiter="\t")
gwas = gwas[['rsid','alt','ref','effect_af','beta','se','pvalue']]
gwas.columns = ['SNP','A1','A2','FRQ','BETA','SE','P']
gwas.to_csv(gwas_file, sep="\t", index = False, compression='gzip')
os.system('mkdir gwas_sim{0}; mv {1} gwas_sim{0}'.format(sim_num, gwas_file))

####### QTL
if (create_predictdb):

	dat = pd.read_csv('{0}/{1}'.format(eqtl_data, eqtl_file), delimiter="\t")

	# make required directories if they don't exist
	if not os.path.isdir('{0}/model_training/dbs'.format(predictdb_pipeline_location)):
		make_predictdb_directory_structure(predictdb_pipeline_location)

	# locuscompare plot to check p-values
	datl = dat[['rsid', 'pvalue']]
	datl.columns = ['SNP', 'eqtl_p']
	gwasl = gwas[['SNP','P']]
	gwasl.columns = ['SNP', 'gwas_p']
	lc = gwasl.merge(datl, on = 'SNP')
	lc.gwas_p = -np.log(lc.gwas_p)
	lc.eqtl_p = -np.log(lc.eqtl_p)
	plt.plot(lc.gwas_p, lc.eqtl_p, 'bo')
	plt.xlabel('-log10 GWAS p-value')
	plt.ylabel('-log10 eQTL p-value')
	plt.title('Simulation {0}'.format(sim_num))
	plt.savefig('sim{0}_locuscompare.png'.format(sim_num))
	plt.clf()

	# create expression file
	expr = pd.read_csv('{0}/{1}'.format(eqtl_data, eqtl_phenotypes), delimiter="\t")
	expr = expr.iloc[:,6:]
	expr.insert(0,"NAME", "ENSG00000000")
	expr.to_csv('tmp.expression{0}.txt'.format(sim_num), sep="\t", index = False)

	# create gene annotation file
	seqname = np.unique(dat.chr)[0]
	gene_id = 'ENSG00000000'
	gene_name = 'nameless_gene'
	start = np.min(dat.snp_pos) + (np.max(dat.snp_pos) - np.min(dat.snp_pos))/2
	end = np.max(dat.snp_pos) + (np.max(dat.snp_pos) - np.min(dat.snp_pos))/2
	gene_type = 'gene'
	with open('tmp.annot{0}.gtf'.format(sim_num), 'w') as annot:
		annot.write('chr\tgene_id\tgene_name\tstart\tend\tgene_type\n')
		annot.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(seqname, gene_id, gene_name, start, end, gene_type))

	# create SNP annotation file
	snp_annot = dat[['chr', 'snp_pos', 'variant_id', 'ref', 'alt', 'effect_af', 'rsid', 'rsid']]
	snp_annot.columns = ['chromosome','pos','varID','ref_vcf','alt_vcf', 'MAF', 'rsid','rsid_dbSNP150']
	snp_annot.insert(5, 'R2', 1)
	snp_annot.to_csv('tmp.snp_annot{0}.chr{1}.txt'.format(sim_num, seqname), sep="\t", index = False)

	# create genotype file (separated by chromosomes)
	os.system('{0}/vcftools --gzvcf {1}/{2} --012 --out {3}'.format(vcftools_location, eqtl_data, eqtl_genotypes, 'tmp.genotypes{0}'.format(sim_num))) 
	## 012 format only outputs chr pos coordinates and not variant id
	## can use snp_annot above to annotate VariantIDs for the next step of the pipeline
	genotypes = pd.read_csv('tmp.genotypes{0}.012'.format(sim_num), delimiter="\t", header = None, index_col = 0)
	indivs = pd.read_csv('tmp.genotypes{0}.012.indv'.format(sim_num), header = None)
	variants = pd.read_csv('tmp.genotypes{0}.012.pos'.format(sim_num), delimiter="\t", header = None)
	# get list of row indices of snp_annot (eqtl summary stats) that match variant list in variants, this is a slightly redundant check due to structured nature of simulated data
	# dat = dat.assign(mergecol = lambda x: ["chr" + str(i) + "_" + str(j) for i,j in zip(dat.chr, dat.snp_pos)])
	# snp_annot = snp_annot.assign(mergecol = lambda x: ["chr" + str(i) + "_" + str(j) for i,j in zip(snp_annot.chromosome, snp_annot.pos)])
	# variants = variants.assign(mergecol = lambda x: [str(i) + "_" + str(j) for i,j in zip(variants[0], variants[1])])
	# snp_annot = snp_annot.set_index('mergecol')
	# snp_annot = snp_annot.reindex(index = variants['mergecol'])
	# snp_annot = snp_annot.reset_index()
	# write genotypes dataframe to genotype file
	genotypes = genotypes.T
	genotypes.insert(0, "varID", list(snp_annot.varID))
	indivs = list(indivs[0])
	indivs.insert(0, "varID")
	genotypes.columns = indivs
	genotypes.to_csv('tmp.genotypes{0}.chr{1}.txt'.format(sim_num, seqname), sep="\t", index = False)
	# os.system('echo -e "\n" >> tmp.genotypes{0}.chr{1}.txt'.format(sim_num, seqname))
	## rsid and R2 are not required, MAF and rsid_dbSNP150 are required

	## Model fitting
	# can set covariates_file to FALSE if we don't want to adjust for covariates
	os.chdir(predictdb_pipeline_location + '/model_training/scripts')
	os.system('Rscript gtex_tiss_chrom_training.R sim{0} {1} {2}/tmp.snp_annot{0}.chr{1}.txt {2}/tmp.annot{0}.gtf {2}/tmp.genotypes{0}.chr{1}.txt {2}/tmp.expression{0}.txt FALSE'.format(sim_num, seqname, working_directory))

	## Create SQL database with weights
	select_signif = False # select only significant weights based on zscore and rho_avg (model fitting parameters)
	os.system('Rscript make_dbs_modified.R sim{0} {1} {2}/tmp.annot{0}.gtf {3} {2}'.format(sim_num, seqname, working_directory, select_signif, working_directory))

## Run PrediXcan
db_name = '../dbs/sim{0}_tw0.5_filtered.db'.format(sim_num)
os.chdir(working_directory)
os.system('{0}/MetaXcan.py --model_db_path {1}/model_training/dbs/sim{2}_tw0.5_filtered.db --covariance {1}/model_training/covariances/sim{2}_nested_cv_chr{3}_covariances.txt --gwas_folder {4}/gwas_sim{2} --gwas_file_pattern ".*gz" --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column P --output_file sim{2}_spredixcan_results.csv'.format(predixcan_location, predictdb_pipeline_location, sim_num, seqname, working_directory))

# clean up temp files
os.system('rm -r gwas_sim{0}'.format(sim_num))
os.system('rm tmp.annot{0}.gtf'.format(sim_num))
os.system('rm tmp.expression{0}.txt'.format(sim_num))
os.system('rm tmp.genotypes{0}.*'.format(sim_num))
os.system('rm tmp.snp_annot{0}.chr*'.format(sim_num))

def make_predictdb_directory_structure(predictdb_pipeline_location):
	os.system('mkdir {0}/model_training/analysis'.format(predictdb_pipeline_location))
	os.system('mkdir {0}/model_training/covariances'.format(predictdb_pipeline_location))
	os.system('mkdir {0}/model_training/dbs'.format(predictdb_pipeline_location))
	os.system('mkdir {0}/model_training/summary'.format(predictdb_pipeline_location))
	os.system('mkdir {0}/model_training/weights'.format(predictdb_pipeline_location))
	os.system('mkdir {0}/prepare_data/covariates'.format(predictdb_pipeline_location))

