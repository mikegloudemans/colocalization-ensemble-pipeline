# Colocalization pipeline
## Check for colocalization between GWAS and eQTL data using eCAVIAR

#### Getting started
First, clone the repository to a location where you can run it on a computing cluster

```
git clone <distribution repository>
```
The pipeline has been tested using Python 2.7.12, and should work for any version of Python 2.7.x. It has not been tested using Python 3.x.

Please ensure you have the following packages, or install them using `pip install`
1. numpy
2. pandas
3. matplotlib

Please check that `/bin` contains the following executables
1. tabix
2. ecaviar
3. vcftools
4. plink

#### Data format in GWAS top hits, full GWAS and full eQTL files
The GWAS top hits, full GWAS and full eQTL files are tab delimited and contain the following fields in sequence
```
chr     rsid    a1      a2      snp_pos beta    se      pvalue
```
where **chr**=chromosome number, **rsid**=RSID identifier, **a1**=Allele 1, **a2**=Allele 2, **snp_pos**=SNP position (build 37), **beta**=Effect size (can be an odds ratio as well, and the corresponding column name should be **or**), **se**=Standard error, **pvalue**=p-value. It is important to keep the order consistent within and across files.

#### Running tabix on full GWAS and eQTL files
Tabix is a tool that greatly increases the speed with which large tab-delimited files can be searched. The pipeline runs on tabix'd files to save execution time. 

Briefly, tabix uses a pair of files to get information - the original compressed file that contains the desired information (compressed using bgzip, not gzip!), and a .tbi index file that tabix uses to rapidly search the first compressed file.

Before running the pipeline, please run tabix on the files that contain the full GWAS and full eQTL data. This needs to be done only once for each large file. Tabix need not be run on the file containing GWAS top hits because it is much smaller. In order to compress each file and create a tabix index, run 

```
./tabix_files.sh <path-to-GWAS-file> <path-to-eQTL-file>
```
This results in compressed and indexed files containing the full GWAS and eQTL information.

#### Directory structure
The main pipeline directory contains 5 subdirectories

1. bin - contains executables required for the pipeline to run
2. data - contains input data including full eQTL/GWAS data and GWAS summary statistics
3. output - contains colocalization results
4. scripts - contains scripts required to run the pipeline
5. tmp - contains temporary files
6. previous_results - general folder for all previous results; old results are moved here automatically 

```
1. bin   
2. data 
  - eqtl_full
  - gwas_full
  - gwas_tophits
3. output 
  - plots
4. scripts                    
5. tmp 
6. previous_results
```

#### Usage
To run the pipeline, ensure all data is placed in the data/ folder. Run

```
source run_colocalization_pipeline.sh <path-to-gwas-tophits-file> <path-to-gwas-file> <path-to-eqtl-file>
```
#### Output subdirectory structure and file descriptions
The file containing posterior colocalization probabilities in all tissues is located in the `output/` directory and has the name *_clpp_status.txt. The format of this tab delimited file is

```
<chr_pos> <tissue_prefix> <gene> <snps.shape> <clpp>
```

where the last column <clpp> is the posterior colocalization probability of GWAS and eQTL signals at the GWAS SNP <chr_pos> for <gene> indicated in the third column in tissue <tissue_prefix>. <snps.shape> is the number of SNPs that were tested for colocalization this site in this tissue.

`output/plots` contains 2 different types of plots at each GWAS top hit. One is a plot of (log eQTL p-value) by (log GWAS p-value), and a colocalization can be observed if certain SNPs have both high eQTL and GWAS p-values. The other is a Manhattan plot of the region around the GWAS top hit, which is indicated by a vertical dashed line. The SNPs are colored by location, and each SNP has the same color on both plots.

