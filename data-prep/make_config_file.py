import json
import glob
import copy

with open("snakemake/colocalization-config.json") as f:
	config = json.load(f)

#trait_set_1 = glob.glob("/oak/stanford/groups/smontgom/mgloud/projects/gwas-download/gwas-download/munge/munged/2021-complete-compendium/hg38/H*/*.gz")
trait_set_1 = glob.glob("/oak/stanford/groups/smontgom/mgloud/projects/gwas-download/gwas-download/munge/munged/2021-complete-compendium/hg38/BMI_Pulit*/*.gz")
trait_set_2 = glob.glob("/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/eqtls/*.gz")
			# + glob.glob("/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/eqtls/*.gz")
example = config["studies"]["Whole_Blood_global"]

config["studies"] = {}

for ts1 in trait_set_1[:1]:
	trait = ts1.split("/")[-1].replace(".txt.gz", "")
	config["studies"][trait] = copy.copy(example)
	config["studies"][trait]["file"] = ts1
	config["studies"][trait]["ref_genome"] = "tg"
	config["studies"][trait]["trait_set"] = 1

for ts2 in trait_set_2[:1]:
	trait = ts2.split("/")[-1].replace(".txt.gz", "")
	config["studies"][trait] = copy.copy(example)
	config["studies"][trait]["file"] = ts2
	config["studies"][trait]["ref_genome"] = "tg"
	config["studies"][trait]["trait_set"] = 2

print(config)

with open("snakemake/whr-test-config.json", "w") as w:
	json.dump(config, w)
