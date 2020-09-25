# Snakemake implementation of colocalization pipeline

## 1. Installation

...


## N. Run the pipeline

Execute the pipeline with:
```bash
conda activate coloc_wrapper
snakemake -j 999 --snakefile Snakefile \
					--cluster-config ds.slurm \
					--latency-wait 90 \
					--configfile config/${tissue}.config \
          --group-components coloc=1000 finemap=1000 \
					--cluster \
					"sbatch --account={cluster.account} \
						--time={cluster.time} \
						--mem={cluster.mem} \
						--cpus-per-task={cluster.nCPUs} \
						--output={cluster.output} \
						--mail-type={cluster.mail}"
```
