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
where **chr**=chromosome number, **rsid**=RSID identifier, **a1**=Allele 1, **a2**=Allele 2, **snp_pos**=SNP position (build 37), **beta**=Effect size (can be an odds ratio as well, and the corresponding column name should be **or**), **se**=Standard error, **pvalue**=p-value.

#### Running tabix on full GWAS eQTL files
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

###### Subdirectory structure
```
1. bin   
2. data 
  - eqtl_full
  - gwas_full
  - gwas_tophits
3. output 
  - coloc_status
  - gwas_by_eqtl_scatterplots
  - snp_overlaps
4. scripts                    
5. tmp   
```
#### Adjustable parameters
One of the main settings that affects colocalization results is the p-value threshold set for eQTL and GWAS data. The pipeline detects colocalization events if for a given SNP, both eQTL and GWAS data p-values exceed their respective thresholds. There might be various reasons to adjust these thresholds (for example, sample size variation between eQTL studies for colocalization analyses and comparison across studies). To adjust GWAS and eQTL p-value thresholds,

#### Usage
To run the pipeline, ensure all data is placed in the data/ folder. Run

```
python find_colocalizations.py <chromosome> <gwas_position> <path-to-gwas-file> <path-to-eqtl-file>
```
#### Output subdirectory structure and file descriptions
1. coloc_status
2. gwas_by_eqtl_scatterplots
3. snp_overlaps

The file containing posterior colocalization probabilities in all tissues is located in the output/ directory and has the name *_clpp_status.txt. The format of this tab delimited file is

```
<gwas_chrom>_<gwas_pos> <tissue_prefix> <gene> <snps.shape> <clpp>
```

where the last column <clpp> is the posterior colocalization probability of GWAS and eQTL signals at the GWAS SNP <gwas_chrom>_<gwas_pos> for <gene> indicated in the third column in tissue <tissue_prefix>. 

 #### Plotting eQTL vs GWAS plots
 output_coloc_plots.py plots log(eQTL p-values) against log(GWAS p-values) for given gene. Run

 ```
python output_coloc_plots.py <chromosome> <gwas_position> <gwas_file> <eqtl_gene>
 ```
where <chromosome> is the chromosome number, <gwas_position> is the position of the top GWAS hit of interest for colocalization, <gwas_file> is the path to the full GWAS data file, and <eqtl_gene> is the ENSEMBL ID (ENSG00000...) of the gene of interest with respect to the current GWAS hit.
