#!/usr/bin/python
# Author: Mike Gloudemans
# Date: 9/27/2017

#
# Run colocalization analysis for all UK Biobank traits.
#

import subprocess
import json
import os

# TODO: Add the option to just analyze one trait, using its pheno-code.

# NOTE: Right now GWAS sample size is missing, and this will only work for methods that don't require it.
# TODO: Make it automatically inject the sample size, case/control ratios, and other metrics needed for COLOC and other methods.
config_template = '''
{
        "eqtl_experiments":	{"/users/mgloud/projects/brain_gwas/data/eqtls/tabix/Whole_Blood_Analysis.cis.eqtl.gz": {"N": 338, "ref": "gtex_full"}},
        "gwas_experiments": {},
	"gwas_threshold": 1e-10,
	"methods": {
		"finemap":{}
	},
	"ref_genomes": {
		"gtex_full": {
			"file": "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/genotypes/OMNI_arrays/variant_calls/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz",
			"af_attribute": "EXP_FREQ_A1",
			"N": 450
		},
		"uk10k": {
			"file": "/users/mgloud/projects/brain_gwas/data/ref/uk10k/chr{0}.vcf.gz",
			"ac_attribute": "AC",
			"N": 1927
		}
	}

}
'''
 
config = json.loads(config_template)
print config

# Inject the GWAS of interest (from the UK Biobank folder) into the 
gwas_list = [g for g in os.listdir("/mnt/lab_data/montgomery/shared/datasets/ukbb/gwas") if g.endswith("assoc.sorted.tsv.gz")]
for g in gwas_list:
    config["gwas_experiments"]["/mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/{0}".format(g)] = {"ref": "uk10k"}

with open("/users/mgloud/projects/brain_gwas/data/config/ukbb.config", "w") as w:
    json.dump(config, w)

# Then run coloc pipeline.
subprocess.check_call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/brain_gwas/data/config/ukbb.config", shell=True)
