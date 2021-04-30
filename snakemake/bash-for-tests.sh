snakemake -j 999 --snakefile skeleton.snakefile \
					--cluster-config test.config \
					--latency-wait 90 \
					--configfile config.json \
					--group-components coloc=500 finemap=500 \
					--cluster \
					"sbatch --account={cluster.account} \
						--partition={cluster.partition} \
						--time={cluster.time} \
						--mem={cluster.mem} \
						--cpus-per-task={cluster.nCPUs} \
						--mail-type={cluster.mail}"


snakemake -j 1 --snakefile skeleton.snakefile --configfile config.json --forceall

snakemake -j 1 --snakefile skeleton.snakefile --configfile config.json --forceall | dot -Tpdf > ~/dag.pdf

snakemake -j 999 --snakefile skeleton.snakefile \
					--cluster-config test.config \
					--latency-wait 90 \
					--configfile config.json \
					--group-components COLOC=500 FINEMAP=500 \
					--cluster \
					"sbatch --account={cluster.account} \
						--partition={cluster.partition} \
						--time={cluster.time} \
						--mem={cluster.mem} \
						--cpus-per-task={cluster.nCPUs} \
						--mail-type={cluster.mail}"
