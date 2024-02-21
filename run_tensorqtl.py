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
print(f'Pandas {pd.__version__}')
print(f'tensorqtl {tensorqtl.__version__}')

outdir='eqtl_nominal'
if not os.path.exists(outdir):
    os.makedirs(outdir)

exprfiles=glob.glob('eqtl/*.expr.bed.gz')
regs = [os.path.basename(file).split('.')[0] for file in exprfiles]
exprdict = {os.path.basename(file).split('.')[0]: file for file in exprfiles}
for reg in regs:
    expres = exprdict[reg]
    covar=expres.replace("expr.bed.gz", "covars.txt")
    if not op.exists(covar):
       raise FileNotFoundError(error_message)("Covars file "+covar+" not found!")

plink =  'genotypes/adhd_gt_n575_maf05'
if op.exists(plink+'.fam'): 
    print("Plink fam file found.")
else:
    print("Plink fam "+plink+".fam not found!")
# ---
pr = genotypeio.PlinkReader(plink)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
print("Genotype dimensions:", end='')
print(genotype_df.shape)
## rename gt columns in case it has 3-digit BrNums
brpat=r'Br(\d{3})$'
gtcols=genotype_df.columns.tolist()
recols=[re.sub(brpat, lambda match: f'Br0{match.group(1)}', s) for s in gtcols]
genotype_df.columns = recols
## Fix chromosomes (add the "chr" prefix) if needed:
if not variant_df.chrom.iloc[0].startswith('chr'):
   variant_df.chrom = [ 'chr' + chrom for chrom in variant_df.chrom]
## select chromosomes - to make sure we have the same chromosomes in our data for each expression dataset
variant_chrom = set(variant_df.chrom)

for reg in regs:
    print(f" Processing region: {reg}")
    tag = reg + '_gwnom' 
    expres = exprdict[reg]
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
