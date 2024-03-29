# Colocalization pipeline

## Author: Mike Gloudemans

### About

This pipeline colocalizes eQTL and GWAS signals using various colocalization methods.
Currently, we have support for COLOC and for the eCAVIAR method, for which we use a
time-efficient implementation, FINEMAP, for the fine-mapping step.

## Table of Contents
- [Setup](#setup)
- [Preparing Input Files](#preparing-input-files)
    - [GWAS and eQTL summary statistics](#gwas-and-eqtl-summary-statistics)
    - [Reference VCF](#reference-vcf)
- [Creating the configuration file](#creating-the-configuration-file)
    - [`eqtl_experiments`](#eqtl_experiments-required)
    - [`gwas_experiments`](#gwas_experiments-required)
    - [`ref_genomes`](#ref_genomes-required)
    - [`methods`](#methods)
- [Running the pipeline](#running-the-pipeline)
- [Pipeline output](#pipeline-output)
    - [FINEMAP/eCAVIAR](#finemap--ecaviar)
    - [COLOC](#coloc)
    - [Errors and skipped variants](#errors-and-skipped-variants)
- [Tips for power users](#tips-for-power-users)
- [Troubleshooting](#troubleshooting)
- [Colocalization software used](#colocalization-software-used)
- [Contributors and acknowledgements](#contributors)

## Setup

To run this pipeline, first clone this repository into your desired directory.

You will need to install [PLINK 1.9](https://www.cog-genomics.org/plink2), [FINEMAP](http://www.christianbenner.com/),
and [tabix](http://www.htslib.org/doc/tabix.html). Binary executables for these tools should be placed in a 
`/bin` directory on the level under the project directory. The PLINK executable must be named `plink`, the FINEMAP executable must be named `finemap`,
and the `tabix` executable must be named `tabix`.

The following python2 modules must also be installed:  
  - [`progress.bar`](https://pypi.org/project/progress/)
  - <add>
  
The following R packages (vN.N or higher) must also be installed:  
  - `coloc` to run COLOC 
  - <add> 

## Preparing Input Files

### GWAS and eQTL summary statistics

GWAS and eQTL input files must be tab-delimited lists of SNPs, containing at least the following columns:

* `chr`: Chromosome
* `snp_pos`: Position of SNP in hg19
* `alt`: Effect allele
* `ref`: Non-effect allele
* `beta`: Estimated effect size of the `alt` allele, when compared with the `ref` allele
* `se`: Standard error of `beta`

QTL files must additionally contain a `gene` or `feature` column that specifies which gene's expression
is being tested. (Multiple features can be given in the same file.)   

QTL files may include the following additional optional columns:  
* `effect_af_eqtl`: Allele frequency of the effect (`alt`) allele  

Formatting GWAS files may be easier using the scripts found at (https://github.com/mikegloudemans/gwas-download).

The column order does not matter, but the column names _must_ be specified in the header of the file.

An example GWAS input file:
```
chr     snp_pos alt     ref     beta        se          pvalue
1       751756  C       T       .013006     .017324     .4528019
1       752566  A       G       -.005243    .0157652    .7394597
1       752721  G       A       -.003032    .0156381    .8462652
1       752894  C       T       .00464      .0162377    .7750657
[...]
```

Depending on which tool(s) you're running, additional columns may be required.
This will be described below.  

GWAS files may include the following additional optional columns:  
* `effect_af_gwas`: Allele frequency of the effect (`alt`) allele  

It is okay to have additional columns beyond the required or optional ones; these will be ignored when the program is run.
What's important is that all the required columns are included.

The GWAS, eQTL, and VCF files must be sorted first by chromosome and then by snp position,
and finally zipped with `bgzip` and indexed using `tabix`. This ensures that the pipeline will run
at a reasonable speed.

The following two commands could be used to index a file named `my_eqtl_file.txt` with the chromosome number
in the second column (corresponds to `-s` parameter) and the chromosome position in the third column (`-b` and `-e` parameters):

```
bgzip -f /home/my_eqtl_file.txt
tabix -f -S 1 -s 2 -b 3 -e 3 /home/my_eqtl_file.txt.gz
```

### Reference VCF

FINEMAP requires a reference VCF to compute the LD between positions of interest.

COLOC requires a reference VCF to estimate MAFs for the variants, unless the MAFs are directly
specified in the GWAS and/or eQTL input files. 

We recommend using the 1000 Genomes VCF, which can be obtained freely from the 1000 Genomes
website in its hg19 or hg38 form.

## Creating the configuration file

The configuration file specifies the specific colocalization methods to be performed,
the parameter settings for these methods, the input and output file locations, and the 
basis for selecting the loci to be tested.

The file must be in JSON format. The top level is a JSON object in which the following parameters may be specified.
To get an idea of what this file should look like, ask for an example config file.

### `out_dir` (recommended)  

The root directory for output files generated by this analysis. 
Default: `/users/mgloud/projects/brain_gwas`  

### `out_dir_group` (recommended)

A directory in which to place all output files generated by this analysis
(relative to the `out_dir` parameter).
If not specified, the output directory will simply be a combination of the config file
name a time stamp of the time at which the pipeline was started running.
 
### `tmp_dir` (recommended)  

The directory where you'll store temporary files (these files will be removed at the end of each run).  
Default: `/users/mgloud/projects/brain_gwas/tmp`  

### `selection_basis`

The process by which the loci to test are selected. Can be `eqtl`, `gwas`, `both`, `snps_from_list`, or `overlap_loci`. It takes the form of a dictionary:
```
    "selection_basis":
    {
        "eqtl": {},
        "gwas": {},
        "both": {},
        "snps_from_list": "/path/to/snps_from_list.txt",
        "overlap_loci": "/path/to/coloc-tests.txt"
     },
     ...
```

Note only **one** of the above methods would actually be included in the config file. For now, just ask me for help when you get to this part.

### `plot_all`  

Plot all colocalization plots from LocusCompare.  
Possible values are `True` and `False`.  

### `plot_none`  

Do not attempt to make colocalization plots from LocusCompare.  
Possible values are `True` and `False`.  
For now, if you are working with GWAS not included in LocusCompare,  
this value **must** be set to `True` to prevent errors.  

---
 
### `eqtl_experiments` (required)

Specifies a JSON object containing additional objects whose keys are the file 
locations of each eQTL file being tested. These files 
must be bgzipped and tabixed for the pipeline to run (see preparation instructions above in the 
"Preparing Input Files" section). 

The eQTL file objects should contain two additional fields:

#### `ref` (required)

A name for the reference VCF to be used for aligning the directions of the eQTL SNPs and in
some cases for estimating the minor allele frequencies. For now, it's best to have the
same reference genome for all eQTL and all GWAS files (hopefully in the future this
will be changed). This name can be any valid string; you will later define the location
and the properties of this reference VCF file in the `ref_genomes` parameter.

If no custom reference VCF is desired, we recommend using the 1K genomes VCF, which
is freely available on the consortium website in an already-indexed form.

#### `eqtl_format` (required)

For now, there are two options available.

Set this parameter to `"effect_size"` if your file contains a column with a signed beta
for each tested SNP and an `se` column specifying the standard error
of the estimated effect.

Set this parameter to `"pval_only"` if your file contains a p-value for each tested SNP   
and one of the following:  
* An `effect_direction` column with `+` and `-` values telling whether the direction of the 
alt allele's effect on gene expression is positive or negative.  
* A `beta` column with signed beta values for each tested SNP.  

---

### `gwas_experiments` (required)

Specifies a JSON object containing additional objects whose keys are the file 
locations of each GWAS file being tested. These files 
must be bgzipped and tabixed for the pipeline to run (see preparation instructions above in the 
"Preparing Input Files" section). 

The GWAS file objects should contain two additional fields:

#### `ref` (required)

A name for the reference VCF to be used for aligning the directions of the eQTL SNPs and in
some cases for estimating the minor allele frequencies. For now, it's best to have the
same reference genome for all eQTL and all GWAS files (hopefully in the future this
will be changed). This name can be any valid string; you will later define the location
and the properties of this reference VCF file in the `ref_genomes` parameter.

If no custom reference VCF is desired, we recommend using the 1K genomes VCF, which
is freely available on the consortium website in an already-indexed form.

#### `gwas_format` (required)

For now, there are two options available.

Set this parameter to `"effect_size"` if your file contains a column with a signed beta or with an
odds ratio for each tested SNP, and an `se` column specifying the standard error
of the estimated effect.

Set this parameter to `"pval_only"` if your file does not contain the above columns, but it contains
a p-value for each tested SNP, and a "direction" column telling whether the direction of the 
alt allele's effect on gene expression is positive or negative.

---

### `ref_genomes` (required)

Specifies a JSON object containing a list of all reference population VCFs used in this analysis.
The reference VCFs are used for calculating LD and sometimes for estimating minor allele frequencies,
and their names should correspond to the reference genomes specified in the `eqtl_experiments` field.

Each reference VCF should be a JSON object with 3 required fields:

#### `file` (required)

The bgzipped, tabixed file where the reference genome is located.

NOTE: Some reference VCFs are given in one long file, while others are given by chromosome. If yours
are given by chromosome, you can include the a `{0}` placeholder in this string denoting the chromosome number
For example, in the filename `1KG_ALL.chr{0}.phase3.vcf.gz`, the colocalization script
will automatically replace the placeholder with the relevant chromosome number when running
the analysis.

#### `af_attribute` (required)

The attribute in the `INFO` field of the VCF that specifies the allele frequency.
This is done so we can filter out low-frequency variants to avoid generating files that
are impossible for colocalization methods to analyze due to edge cases.

For example, in 1K genomes VCFs the info field contains the following string:
```
AC=15;AF=0.00299521;AN=5008;NS=2504;DP=2036;EAS_AF=0.002;AMR_AF=0;AFR_AF=0.0091 ...
```

We want to extract the allele frequency, which is the second parameter listed here, so we
enter `"AF"` for the `af_attribute` field.

(For now, if you're having trouble figuring this out, just ask me. This part isn't very well-polished
right now.)

#### `N` (required)

The total number of individuals profiled in the reference population.
(For the entire 1K genomes VCF, this number will be 2504.)

---

### `methods`

An object containing a list of the colocalization methods to use for this run.
Each entry is itself a JSON object whose key is the method to run and whose value
is an object specifying the parameter setting for this colocalization method.

Right now there are only two tested colocalization methods available: `coloc`
and `finemap`.

Example: to run both FINEMAP and COLOC with no parameter modifications, we would specify
```
{
    ...,
    "methods":
    {
        "finemap":{},
        "coloc":{}
    },
    ...
}
```

We now describe additonal parameters and setting needed to run
each colocalization method:

---

### `coloc`

The COLOC pipeline has no modifiable parameters; however, in
order to run this method, some additional properties need
to be specified.

Each GWAS in the `gwas_experiments` section must specify not only a `ref` and a
`gwas_format`, but also the following parameters.

#### `type`

Either `"cc"` or `"quant"`, depending on whether the GWAS is a quantitative or case-control GWAS.  
If the same `type` applies to all of the files in `gwas_experiments`, you may instead ignore  
the `gwas_experiments`-level `type` value and instead assign `"cc"` or `"quant"` to the  
global parameter `"gwas_type"`. 

#### `gwas_experiments: N`

The sample size of the GWAS. Currently, I'm not sure whether it's valid to estimate the
sample size if it's not known; I hope to investigate this eventually.

Additionally, to run COLOC, the `eqtl_experiments` need to have this parameter added:

#### `eqtl_experiments: N`

The total sample size of the eQTL study. 

Alternatively, if you are performing colocalizations for a large number of GWAS and QTL files  
and do not want to specify `"N"` within the config file, you may define the following global  
parameters:  

* `gwas_sample_sizes`  
* `eqtl_sample_sizes`  

Both of these parameters may to tab-delimited files with no header and two columns:  
first, the file name (not the full path); second, the sample size of the corresponding study.  
For example, `gwas_sample_sizes` may point to a file that looks like this:  
```
coloc_imputed_ADIPOGen_Adiponectin.txt.gz	29304
coloc_imputed_Astle_et_al_2016_Eosinophil_counts.txt.gz	173480
coloc_imputed_Astle_et_al_2016_Granulocyte_count.txt.gz	173480
coloc_imputed_Astle_et_al_2016_High_light_scatter_reticulocyte_count.txt.gz	173480
coloc_imputed_Astle_et_al_2016_Lymphocyte_counts.txt.gz	173480
```
The file names must exactly match the corresponding `basename` of the files specified in   
`gwas_experiments` and `eqtl_experiments`.  

Look [here](https://github.com/mikegloudemans/brain_gwas/blob/coloc_dev/config/coloc_example.config) for an example of a config file written for COLOC.  

----

### `finemap`

Currently, the FINEMAP pipeline has no modifiable parameters.


## Running the pipeline

Once you have specified your desired options in the configuration file, it's
simple to run the pipeline.

```
python ./dispatch.py [config_file] [num_threads]
```

If you are running the pipeline from SCG, load the following modules first:
```bash
module load r/3.6 # R v3.6
module load miniconda/2 # python2
module load plink/1.90b6.13 # PLINK v1.90b
module load tabix # tabix
```

## Pipeline output

The results of the pipeline are placed in a time-stamped folder in the output directory.

[LocusCompareR](https://github.com/boxiangliu/locuscomparer) visualizations of significant loci are 
available in the `plots` subdirectory. Please see the LocusCompareR website for more details on
how to interpret these plots if you are unfamiliar with LocusZoom-style plots.

The rest of the output will vary depending on which colocalization tools you ran. The outputs
for specific tools are described below.

### FINEMAP / eCAVIAR

At the top level, the directory contains a file `[gwas_name]_finemap_clpp_status.txt`.
This file contains a row for each locus tested. The sixth column is the CLPP score.
For more information on interpreting the CLPP score, see [Hormozdiari et al](https://www.ncbi.nlm.nih.gov/pubmed/27866706).
We suggest filtering tested loci by the CLPP score and then visualizing the locus for further confirmation, to avoid
false positives.

A good starting filter is to select loci with at least 100 variants tested, a CLPP score of > 0.02,
and a top -log(pvalue) of > 5 for both the eQTL and the GWAS studies. These filters can be made more
stringent as needed if too many false positives are being detected.

The `finemap` directory contains two file for each colocalization test, showing the posterior causal
probabilities for each SNP at the locus. These probabilities are shown in the `snp_prob` column.

### COLOC

A file will be created in the top directory with the form `[gwas_name]_coloc_h4pp_status.txt`.
There will be an `h4pp` column in this file that contains the posterior probability of 
H4 (colocalization) as described the COLOC model framework. 

### Errors and skipped variants

Regardless of which colocalization methods were selected, two files will be created
named `ERROR_variants.txt` and `skipped_variants.txt` that aid in troubleshooting
failed colocalization runs.

`ERROR_variants.txt` contains a list of loci for which the pipeline crashed before
running to completion. A Python trackback report will be listed in the last column
whenever possible, and will hopefully be helpful in diagnosing the problem; if it
is unclear though, you're welcome to contact me or to file an issue report on the
Github page.

`skipped_variants.txt` contains variants that were skipped for known reasons rather
than due to unhandled issues in the pipeline. Some common reasons for skipping variants
are when no variants overlap between the GWAS, eQTL, and VCF files, or when the width
of the overlapping region is too small to perform a proper colocalization analysis.
In some rare cases we also find variants that are skipped because no genes were tested
for the eQTL study in that region; when this happens the description is listed as
"Gene desert".

## Tips for power users

### Assembling complex configuration files

In some cases, if you want to test a very large number of GWAS files (e.g. UK Biobank),
it can be cumbersome to actually create the config file required to analyze all the
files. Similar, in some cases you might only want to test certain combinations of
GWAS and eQTL files at certain loci, rather than every possible pairwise combination.

If this is true for you, I recommend that instead of manually assembling the config file,
you load in a template of the JSON file and construct it programmatically using Python's
`json` library or your language of choice.

I hope to eventually make available some code I've written for performing this in
a principled and automated manner on any arbitrary set of GWAS and eQTL files according
to the combinations passing certain thresholds.

## Troubleshooting

My goal is to make these pipelines as easily runnable and as intuitive as possible!
Please contact me with any problems you encounter, or post an issue report, and
I'll happily work with you to help you out or to fix any bugs you might find.

## Colocalization software used

[FINEMAP v1.1](http://www.christianbenner.com/)
cite

[add others before making public]

### Contributors

* Mike Gloudemans
* Abhiram Rao
* Brunilda Balliu
* Nicole Gay
* Boxiang Liu
* Jeremy Tien


### Additional Acknowledgments

A big thanks to the following people with help troubleshooting and testing the pipeline:

* Pagé Goddard
* Derek Klarin
