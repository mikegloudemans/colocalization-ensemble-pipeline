# cluster execution
snakemake -j 100 --cluster-config slurm-config.json \
					--latency-wait 90 \
					--cluster \
					"sbatch --account={cluster.account} \
						--partition={cluster.partition} \
						--time={cluster.time} \
						--mem={cluster.mem} \
						--cpus-per-task={cluster.nCPUs} \
						--mail-type={cluster.mail} \
						--output={cluster.output}" 

# dry run
snakemake -n --configfile colocalization-config.json

# execute on current node
snakemake -j 1 --configfile colocalization-config.json

# print DAG to home
# this won't work if the snakefile prints anything to stdout
snakemake -j 1 --configfile colocalization-config.json --forceall | dot -Tpdf > ~/dag.pdf

