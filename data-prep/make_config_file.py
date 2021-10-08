import json
import glob
import copy

with open("snakemake/colocalization-config.json") as f:
	config = json.load(f)

'''
ts1_unglobbed = ["/oak/stanford/groups/smontgom/mgloud/projects/insulin_resistance/data/gwas/formatted/sumstats/hg38/*/*.gz"]
ts2_unglobbed = ["/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/eqtls/Adipose*.gz",
			"/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/eqtls/Pancreas*.gz",
			"/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/eqtls/Liver*.gz",
			"/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/eqtls/Skeletal*.gz",
			"/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/sqtls/Adipose*.gz",
			"/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/sqtls/Pancreas*.gz",
			"/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/sqtls/Liver*.gz",
			"/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/sqtls/Skeletal*.gz"]
'''
ts1_unglobbed = ["/oak/stanford/groups/smontgom/mgloud/projects/insulin_resistance/data/gwas/formatted/sumstats/hg38/BMI_GIANT_2018/BMI_GIANT_2018.txt.gz"]
ts2_unglobbed = ["/oak/stanford/groups/smontgom/mgloud/gtex-data-indexed/eqtls/Liver*.gz"]

trait_set_1 = []
trait_set_2 = []

for ts1 in ts1_unglobbed:
	trait_set_1 += glob.glob(ts1)

for ts2 in ts2_unglobbed:
	trait_set_2 += glob.glob(ts2)
 
example = config["studies"]["Whole_Blood_global"]

config["studies"] = {}

#for ts1 in trait_set_1[:1]:
for ts1 in trait_set_1:
	trait = ts1.split("/")[-1].replace(".txt.gz", "")
	config["studies"][trait] = copy.copy(example)
	config["studies"][trait]["file"] = ts1
	config["studies"][trait]["ref_genome"] = "tg"
	config["studies"][trait]["trait_set"] = 1
	if "N" in config["studies"][trait]:
		del config["studies"][trait]["N"]
	if "trait_type" in config["studies"][trait]:
		del config["studies"][trait]["trait_type"]

#for ts2 in trait_set_2[:1]:
for ts2 in trait_set_2:
	trait = ts2.split("/")[-1].replace(".txt.gz", "")
	config["studies"][trait] = copy.copy(example)
	config["studies"][trait]["file"] = ts2
	config["studies"][trait]["ref_genome"] = "tg"
	config["studies"][trait]["trait_set"] = 2
	if "N" in config["studies"][trait]:
		del config["studies"][trait]["N"]
	if "trait_type" in config["studies"][trait]:
		del config["studies"][trait]["trait_type"]

print(config)

with open("snakemake/whr-test-config.json", "w") as w:
	json.dump(config, w)
