# Colocalization pipeline
## Check for colocalization between GWAS and eQTL data using eCAVIAR

#### Getting started
First, clone the repository to a location where you can run it on a computing cluster

```
git clone <distribution repository>
```

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
  - eqtl
  - gwas_summary
  - gwas
3. output 
  - coloc_status
  - gwas_by_eqtl_scatterplots
  - snp_pverlaps
4. scripts                    
5. tmp   
```
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
<gwas_chrom>_<gwas_pos> <tissue_prefix> <gene> <conditional_level> <snps.shape> <clpp>
```

where the last column <clpp> is the posterior colocalization probability of GWAS and eQTL signals at the GWAS SNP <gwas_chrom>_<gwas_pos> for <gene> indicated in the third column in tissue <tissue_prefix>. The <conditional_level> is a binary variable that is 0 if the colocalization does not pass the preset threshold, or 1 if it does.

 #### Plotting eQTL vs GWAS plots
 output_coloc_plots.py plots log(eQTL p-values) against log(GWAS p-values) for given gene. Run

 ```
./output_coloc_plots.py <chromosome> <gwas_position> <gwas_file> <eqtl_gene>
 ```
where <chromosome> is the chromosome number, <gwas_position> is the position of the top GWAS hit of interest for colocalization, <gwas_file> is the path to the full GWAS data file, and <eqtl_gene> is the ENSEMBL ID (ENSG00000...) of the gene of interest with respect to the current GWAS hit.
