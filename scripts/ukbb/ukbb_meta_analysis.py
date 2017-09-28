#!/usr/bin/python
# Author: Mike Gloudemans
# Date: 9/27/2017

#
# Run colocalization analysis for all UK Biobank traits.
#

# Let's start with just one trait though.

# NOTE: Right now GWAS sample size is wrong, and this will only be valid for FINEMAP-like methods.
# TODO: Make it automatically inject the sample size, case/control ratios, and other metrics needed for COLOC and other methods.
config_template = '''
{{
        "eqtl_experiments":	{{"/users/mgloud/projects/brain_gwas/data/eqtls/tabix/Whole_Blood_Analysis.cis.eqtl.gz": {{"N": 338, "ref": "gtex_full"}}}},
	"gwas_experiments": {{"{0}": {{"N":150064, "ref", "uk10k"}}}},
	"gwas_threshold": 5e-8,
	"methods": {{
		"finemap":{{}}
	}},
	"ref_genomes": {{
		"gtex_full": {{
			"file": "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/genotypes/OMNI_arrays/variant_calls/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz",
			"af_attribute": "EXP_FREQ_A1",
			"N": 450
		}},
		"uk10k": {{
			"file": "/users/mgloud/projects/brain_gwas/data/ref/uk10k/chr{{0}}.vcf.gz",
			"ac_attribute": "AC",
			"N": 1927
		}}
	}}

}}
'''

gwas = "/users/mgloud/projects/gwas/scripts/data/T14.assoc.sorted.tsv.gz"

print config_template.format(gwas)

with open("/users/mgloud/projects/brain_gwas/data/config/ukbb.config", "w") as w:
    w.write(config_template.format(gwas))

# Then run coloc pipeline.
subprocess.check_call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/brain_gwas/data/config/ukbb.config", shell=True)
