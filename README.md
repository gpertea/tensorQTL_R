R based workflow to prepare RangeSummarizedExperiment data for tensorQTL.

## Workflow

This workflow is designed to prepare RangeSummarizedExperiment data and genotype data for use with the tensorQTL package. There are 4 main steps to this workflow, involving 4 scripts or notebooks. These 4 scripts have to be edited manually to specify the input and output directories, or file name prefixes.

To prepare for this workflow, fetch the scripts from the github repository: git clone <https://github.com/gpertea/tensorQTL_R>

This will create a `tensorQTL_R` directory with the scripts in it. Change to this directory, and for convenience, create (or symlink) adirectory or directories with your input files (RSEs and genotype data) under this directory.

### 1. Prepare genotype data and SNP PCs: 01_calc_snpPCs.R

This script requires `plink` and `plink2` binaries to be available in the PATH of the current shell environment. The script prepares the genotype data by applying some minimal QC and then it calculates the population stratification components ("SNP PCs"). This is done using the `01_calc_snpPCs.R` script. Edit this script to specify the appropriate values for the following variables at the top of the script:

-   **`gtinput`** : the path to the `.vcf.gz` or the plink `.bed` file
-   **`bedout`** : the plink .bed file prefix for the QC-ed genotype data created by the script

The **`bedout`** .bed file prefix will be the one to use as "genotype input" for the `plink` variable in the tensorQTL python notebook (step 3 below, `3_run_tensorqtl.ipynb`).

Note that besides the QC-ed `bedout` plink genotype data file, the other important output file which is needed for step 2 below is the `snpPCs` file. This will be created as ***`bedout`*****`.snpPCs.tab`**

### 2. Prepare expression data for tensorQTL: 02_prep_tensorQTL_RSE_input.R

The second step is to prepare the expression data for tensorQTL. This is done using the `02_prep_tensorQTL_RSE_input.R` script. Edit this script to specify the appropriate values for the following variables at the top of the script:

-   **`fgrse`** : the path to the gene RSE in the variable at the top of the script. assuming this RSE file is produced by SPEAQeasy, the rse_tx, rse_exon and rse_jx files for the other genomic features will be found in the same directory and expression data will be generated for all the features found.
-   **`fsnp_pcs`** : the path to the SNP PCs file created in step 1 above
-   **`modelstr`** : the design model string to use for data preparation, with the covariates to be used. It should be a string that can be used as an argument to the `formula` function in R. For example, if the model is to include the top SNP PCs as covariates (usually the case), the string should be something like `"~Age+snpPC1+snpPC2+snpPC3+snpPC4+snpPC5"`.
-   **`dsname`** : the name of the dataset to be used for these files prepared for tensorQTL. This will be used as the prefix for the output files, which are going to be *`dsname`*`.feature.expr.bed.gz` and *`dsname`*`.feature.covars.txt`. Note that the same `dsname` will be used as the `dataset` variable in the tensorQTL python notebook (step 3 below, `3_run_tensorqtl.ipynb`).
-   **`fgeno2rna`** : if not NULL, this file path should point to a file that maps the genotype IDs (used in the genotypes and the snpPCs file generated at step 1 above) to the sample IDs in the RNA-seq data (the colnames of the RSEs, or SAMPLE_ID column). This is only needed when the genotype IDs in the genotype data are different from the sample IDs in the RNAseq data. The file should be a tab-delimited file with two columns, without a column header line, the first column being the genotype sample IDs and the 2nd column being the sample IDs in the RNAseq data (RSE column names).

### 3. Run tensorQTL: 03_run_tensorqtl.ipynb

In order to start this notebook, you have to activate the proper `tensorqtl` environment:

`conda activate tensorqtl`.

Then jupyter lab can be started with `jupyter lab` command, and chunks have to be executed in order. 
Edit the 2nd cell of the notebook before executing it, to specify the appropriate values for the variables there. 


### 4. Gathering and filtering tensorQTL output: 04_gather_eqtl_results.Rmd

This script is an Rmarkdown notebook that takes the tensorQTL `.parquet` files (generated one per chromosome), calculates FDR and filters to keep only mappings that have FDR\<0.05, and merges the remaining mappings into a single `data.table` object which is saved as a `.qs` file which can later loaded with `qs::qread()` for further analysis.
