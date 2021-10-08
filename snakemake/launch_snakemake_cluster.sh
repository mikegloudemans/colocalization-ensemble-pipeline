# cluster execution
snakemake -j 24 --cluster-config slurm-config.json \
					--latency-wait 90 \
                                        --stats stats.txt \
                                        --keep-going \
					--rerun-incomplete \
					--configfile rna-config.json \
					--cluster \
					"sbatch --account={cluster.account} \
						--partition={cluster.partition} \
						--time={cluster.time} \
						--mem={cluster.mem} \
						--cpus-per-task={cluster.nCPUs} \
						--output={cluster.output}" 
