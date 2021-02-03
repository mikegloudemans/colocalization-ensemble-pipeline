# Snakemake implementation of colocalization pipeline

## 1. Installation

### 1.1 Install `conda` environment  
Install `conda` (python3) if it is not already installed. The [miniconda installer](https://docs.conda.io/en/latest/miniconda.html) is a convenient way to do this.  

Create a new `conda` environment with the correct dependencies. This can take a while, so it is recommended to run this command in a `screen`/`tmux` session. It also requires more than 16G of RAM:  
```bash
conda env create -f environment.yml
```

### 1.2 Install `FINEMAP v1.4`
Follow the instructions [here](http://www.christianbenner.com/).  

### 1.3 Install...(other dependencies)
#TODO


## N. Run the pipeline

Execute the pipeline with:
```bash
conda activate colocalization-ensemble
snakemake -j 999 --snakefile Snakefile \
					--cluster-config slurm.config \
					--latency-wait 90 \
					--configfile config.json \
					--group-components coloc=500 finemap=500 \
					--cluster \
					"sbatch --account={cluster.account} \
						--partition={cluster.partition} \
						--time={cluster.time} \
						--mem={cluster.mem} \
						--cpus-per-task={cluster.nCPUs} \
						--output={cluster.output} \
						--mail-type={cluster.mail}"
```