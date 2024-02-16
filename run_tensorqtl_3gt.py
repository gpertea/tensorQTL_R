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

WIN_SIZE=1000000
MAF_MIN=0.01 # MAF threshold for full genotypes

outdir='eqtl_out'

if not os.path.exists(outdir):
    os.makedirs(outdir)
genotypes=glob.glob('genotypes/adhd_gt_n458*.bed')
exprfiles=glob.glob('eqtl_inputs/*.expr.bed.gz')
## initial pass to check for expr-covar pairings
for expres in exprfiles:
    reg=op.basename(expres).split('.')[0]
    feature=op.basename(expres).split('.')[1]
    covar=expres.replace("expr.bed.gz", "covars.txt")
    # check if the covariate file exists
    if not op.exists(covar):
        print(f" Covariate file {covar} not found!")
        sys.exit(1)

for plink in genotypes:
  if op.exists(plink.replace('.bed', '.fam')): 
      print("Reading plink genotype: "+plink.replace('.bed', ''))
  else:
      print("Plink fam "+plink+".fam not found!")
  gtname='genomewide'
  if (re.search('27', plink)):
    gtname='27_loci'
  else:
     if (re.search('cred', plink)):
       gtname='credvars'
  print("--------- Genotype: "+gtname)
  #continue
  pr = genotypeio.PlinkReader(plink)
  genotype_df = pr.load_genotypes()
  variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
  print("  Genotype data size:", end='')
  print(genotype_df.shape)
  ## rename gt columns in case it has 3-digit BrNums
  brpat=r'Br(\d{3})$'
  gtcols=genotype_df.columns.tolist()
  recols=[re.sub(brpat, lambda match: f'Br0{match.group(1)}', s) for s in gtcols]
  genotype_df.columns = recols
  ## Fix chromosomes (add the "chr" prefix) if needed:
  if not variant_df.chrom.iloc[0].startswith('chr'):
     #print("  ..adding 'chr' prefix to chromosome names")
     variant_df.chrom = [ 'chr' + chrom for chrom in variant_df.chrom]
  ## select chromosomes - to make sure we have the same chromosomes in our data for each expression dataset
  variant_chrom = set(variant_df.chrom)
  for expres in exprfiles:
    reg=op.basename(expres).split('.')[0]
    feature=op.basename(expres).split('.')[1]
    covar=expres.replace("expr.bed.gz", "covars.txt")
    tag = reg + '.'+feature + '-' + gtname
    print(f"running tensorQTL on {tag}")
    covariates_df = pd.read_csv(covar, sep='\t', index_col=0).T
    #print("Covariates dim:", end='')
    #print(covariates_df.shape)
    #print(covariates_df.iloc[:5, :5])
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expres)
    #print("Phenotype dimensions:", end='')
    #print(phenotype_df.shape)
    ## use the same chromosome set
    express_chrom = set(phenotype_pos_df.chr)
    ## make sure we have the same chromosomes in our data
    ## print an error message and exit if there is no intersection
    if len(variant_chrom.intersection(express_chrom))==0:
        print(f"  No common chromosomes between genotype and expression data for {tag}")
        sys.exit(1)
        continue

    #assert len(variant_chrom.intersection(express_chrom))>0
    if express_chrom - variant_chrom:
        chrom_filter = phenotype_pos_df.chr.isin(variant_chrom)
        if (len(chrom_filter)<phenotype_df.shape[0]):
          phenotype_df = phenotype_df[chrom_filter]
          phenotype_pos_df = phenotype_pos_df[chrom_filter]
    ## make sure we keep only the genotypes for the expression samples
    cols=phenotype_df.columns.tolist()
    geno_df=genotype_df.loc[:, cols]
    maf_min=0
    ## if genotype_df has more than 600,000 variants, use a MAF threshold of 0.01
    if (genotype_df.shape[0]>600000):
       maf_min=MAF_MIN
       print(f"  MAF threshold set to {maf_min}.")
    else:
       print("  No MAF threshold.")
    ## run tensorQTL:
    cis.map_nominal(geno_df, variant_df, phenotype_df, phenotype_pos_df, prefix = tag, 
           covariates_df= covariates_df, maf_threshold=maf_min, window=WIN_SIZE, output_dir=outdir, 
           verbose=False)
print("\nAll done.")

