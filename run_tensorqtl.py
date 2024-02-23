#!/usr/bin/env python
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis
import time
import os
import os.path as op
import sys
import glob
import re
print(f'PyTorch {torch.__version__}')
#print(torch.__version__) 
print('CUDA available: {} ({})'.format(torch.cuda.is_available(), torch.cuda.get_device_name(torch.cuda.current_device())))
print(f'Pandas {pd.__version__}')
print(f'tensorqtl {tensorqtl.__version__}')

outdir='eqtl_out'
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
in_dir='eqtl_inputs' ## directory path where the tensorqtl input files can be found:
                     ## dsname.gene.expr.bed.gz and dsname.covars.txt
dsname='mdd_amyg' ## dataset name, this is the prefix for input and output file names
#plink='genotypes/gt_bsp12_n114_rIDs' ## plink bed prefix for genotype data
plink='genotypes/mdd_maf01'
mapping_file=None ## comment the line below if genotype IDs are the same with RNAseq sample IDs.
#mapping_file='genoID2rnaID.tab' # genotype IDs will be changed to their RNAseq mappings in the 2nd column
mapping_file='genotypes/mdd_geno2rna.tab'
## if not None, mapping_file must be the same with the one used for 02_prep_tensorQTL_RSE_input.R
col_map = None
if mapping_file is not None and op.exists(mapping_file):
    print(f"mapping genotype IDs to RNAseq sample IDs based on mapping file: {mapping_file}")
    mapping_df = pd.read_csv(mapping_file, sep="\t", header=None, names=["gt_id", "r_id"])
    col_map = dict(zip(mapping_df['gt_id'], mapping_df['r_id']))
if op.exists(plink+'.fam'): 
    print("Plink fam file found.")
else:
    raise Exception("Plink fam "+plink+".fam not found!")

exprfiles=glob.glob(in_dir+'/'+dsname+'.*.expr.bed.gz')
features = [os.path.basename(file).split('.')[1] for file in exprfiles]
exprdict = {os.path.basename(file).split('.')[1]: file for file in exprfiles}
for feat in features:
    expres = exprdict[feat]
    covar=expres.replace("expr.bed.gz", "covars.txt")
    if not op.exists(covar):
       raise FileNotFoundError(error_message)("Covars file "+covar+" not found!")
print("Features found: ",features)
print(" Change the `features` array below to set the features to be processed:")
## CHANGE here and uncomment, if needed, in order to select the features to process:
features = ['gene']
print("Features to process: ", features)

pr = genotypeio.PlinkReader(plink)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
print("Genotype dimensions:", end='')
print(genotype_df.shape)
# if a mapping file is given
# Check if all genotype_df column names are in the mapping file
#all_columns_found = genotype_df.columns.isin(mapping_df["current_column"]).all()

if col_map is not None:
  # Rename columns using the dictionary
  genotype_df.rename(columns=col_map, inplace=True)
## Now genotype_df contains the updated column names
#print(genotype_df.iloc[:5, :7])  # Display the first few rows

##  etc. Fix chromosomes (add the "chr" prefix) if needed:
if not variant_df.chrom.iloc[0].startswith('chr'):
   variant_df.chrom = [ 'chr' + chrom for chrom in variant_df.chrom]
## select chromosomes - to make sure we have the same chromosomes in our data for each expression dataset
variant_chrom = set(variant_df.chrom)

for feat in features:
    print(f" Processing feature: {feat}")
    tag = dsname+'.'+feat
    expres = exprdict[feat]
    covar=expres.replace("expr.bed.gz", "covars.txt")
    covariates_df = pd.read_csv(covar, sep='\t', index_col=0).T
    print("Covariates dim:", end='')
    print(covariates_df.shape)
    #print(covariates_df.iloc[:5, :5])
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expres)
    print("Phenotype dimensions:", end='')
    print(phenotype_df.shape)
    ## use the same chromosome set
    express_chrom = set(phenotype_pos_df.chr)
    assert len(variant_chrom.intersection(express_chrom))>0
    if express_chrom - variant_chrom:
        chrom_filter = phenotype_pos_df.chr.isin(variant_chrom)
        if (len(chrom_filter)<phenotype_df.shape[0]):
          phenotype_df = phenotype_df[chrom_filter]
          phenotype_pos_df = phenotype_pos_df[chrom_filter]
    ## make sure we keep only the genotypes for the expression samples
    cols=phenotype_df.columns.tolist()
    geno_df=genotype_df.loc[:, cols]
    ## run tensorQTL:
    cis.map_nominal(geno_df, variant_df, phenotype_df, phenotype_pos_df, prefix = tag, covariates_df= covariates_df,
                maf_threshold=0.05, window=500000, output_dir= outdir, verbose=False)
print("All done.")
