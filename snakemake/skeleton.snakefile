#!/bin/python3

import sys
import os
import glob

# change working directory
os.chdir(config['outdir'])
# make outdir and tmpdir 
subprocess.call('mkdir -p log/cluster', shell=True)
subprocess.call('mkdir -p {}'.format(config['tmpdir']), shell=True)

# define colocalization methods
methods = [k for k in config['colocalization_methods']]

def get_trait_set(set_number):
	trait_set = set()
	for trait in config['studies']:
		if 'trait_set' not in config['studies'][trait]:
			continue
		if config['studies'][trait]['trait_set'] == set_number:
			trait_set.add(trait)
	return(trait_set)

# assume there are two sets of traits for now
# all studies in trait_set_1 are tested against all studies in trait_set_2
# we should probably also add the option to provide a two-column input of all desired trait comparisons
trait_set_1 = get_trait_set(1)
trait_set_2 = get_trait_set(2)
# trait_set_1 = ['imputed_ICBP_DiastolicPressure']
# trait_set_2 = ['Whole_Blood_local']
studies = list(chain(*[list(trait_set_1), list(trait_set_2)]))

# don't let wildcards contain directory
wildcard_constraints:
    study="[^/]*",
    trait1="[^/]*",
    trait2="[^/]*",
   	method="[^/]*"

# print(studies)
# print(trait_set_1)
# print(trait_set_2)

# a = expand('results/{trait1}.{trait2}.{method}.results.txt', 
# 	trait1 = trait_set_1, 
# 	trait2 = trait_set_2,
# 	method = methods),
# b= expand('preprocess/{study}.formatted.summary_stats.txt.gz',
# 	study = studies)
# print(a)
# print(b)

# quit()

# make dictionary of all traits to file paths 
trait_to_raw_file = {}
for trait in studies:
	trait_to_raw_file[trait] = config['studies'][trait]['file']

if len(trait_set_1) == 0 or len(trait_set_2) == 0:
	sys.exit('Two sets of traits were not defined in the config file. (more detail...)')

# rules that can be run on the login node
localrules: all,stage,merge_results


# here we define colocalization methods and trait combinations 
rule all:
	input:
		expand('results/{trait1}.{trait2}.{method}.results.txt', 
			trait1 = trait_set_1, 
			trait2 = trait_set_2,
			method = methods),
		expand('preprocess/{study}.formatted.summary_stats.txt.gz',
			study = studies)


rule stage:
	input:
		lambda wildcards: trait_to_raw_file[wildcards.study]
	output:
		'summary_stats/{study}.raw.summary_stats.txt.gz'
	shell:
		'''
		ln -s {input} {output}
		'''


rule preprocess:
	input:
		'summary_stats/{study}.raw.summary_stats.txt.gz'
	output:
		'preprocess/{study}.formatted.summary_stats.txt.gz'
	shell:
		'''
		touch {output}
		'''


rule overlap:
	input:
		trait1 = 'preprocess/{trait1}.formatted.summary_stats.txt.gz',
		trait2 = 'preprocess/{trait2}.formatted.summary_stats.txt.gz'
	output:
		'overlap/{trait1}.{trait2}.overlap.txt'
	shell:
		'''
		echo 'iter' > {output}
		rand=$((1 + $RANDOM % 10)) # random number 1-10
		for i in $(seq 1 $rand); do 
			echo ${{i}} >> {output}
		done
 		'''


checkpoint split_loci:
	input:
		'overlap/{trait1}.{trait2}.overlap.txt'
	output:
		directory('split_loci/{trait1}.{trait2}')
	run:
		os.mkdir(output[0])
		# read in header first
		with open(input[0]) as f:
			header = f.readline()
		# write out one file per line
		n=0
		with open(input[0], 'r') as f:
			next(f) # skip header
			for line in f:
				with open('split_loci/{}.{}/chunk{}.txt'.format(wildcards.trait1, wildcards.trait2, n), 'w') as o:
					o.write(header)
					o.write(line)
				n+=1


rule coloc:
	input:
		'split_loci/{trait1}.{trait2}/chunk{n}.txt'
	output:
		'coloc/{trait1}.{trait2}.{n}.results.txt'
	group:
		'COLOC'
	shell:
		'''
		echo {wildcards.trait1} {wildcards.trait2} {wildcards.n} > {output}
		'''


rule finemap:
	input:
		'split_loci/{trait1}.{trait2}/chunk{n}.txt'
	output:
		# make these temporary
		'finemap/{trait1}.{trait2}.{n}.results.txt'
	group:
		'FINEMAP'
	shell:
		'''
		echo {wildcards.trait1} {wildcards.trait2} {wildcards.n} > {output}
		'''


def collect_colocalization_results(wildcards):
	checkpoint_output = checkpoints.split_loci.get(**wildcards).output[0]
	# split_loci/{trait1}.{trait2}

	mfiles = expand("{method}/{trait1}.{trait2}.{n}.results.txt",
		method=wildcards.method,
		trait1=wildcards.trait1,
		trait2=wildcards.trait2,
		n=glob_wildcards(os.path.join(checkpoint_output, "chunk{n}.txt")).n)

	print(mfiles)

	return mfiles


rule merge_results:
	input:
		collect_colocalization_results
	output:
		'results/{trait1}.{trait2}.{method}.results.txt'
	shell:
		'''
		cat {input} > {output}
		# cat {input} > results/{wildcards.trait1}.{wildcards.trait2}.{wildcards.method}.results.tmp.txt
		# # remove all but first header
		# head -1 results/{wildcards.trait1}.{wildcards.trait2}.{wildcards.method}.results.tmp.txt > results/{wildcards.trait1}.{wildcards.trait2}.{wildcards.method}.results.txt
		# # get first word in header
		# first_word=$(sed 's/\s.*$//' results/{wildcards.trait1}.{wildcards.trait2}.{wildcards.method}.results.txt)
		# grep -v "^${{first_word}}	" results/{wildcards.trait1}.{wildcards.trait2}.{wildcards.method}.results.tmp.txt >> results/{wildcards.trait1}.{wildcards.trait2}.{wildcards.method}.results.txt
		# rm results/{wildcards.trait1}.{wildcards.trait2}.{wildcards.method}.results.tmp.txt
		# gzip results/{wildcards.trait1}.{wildcards.trait2}.{wildcards.method}.results.txt
		'''

