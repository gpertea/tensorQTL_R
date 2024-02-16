#!/usr/bin/env python
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis
import time
import sys
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

reg = sys.argv[1]
if not reg or (reg != 'Amygdala' and reg != 'sACC'):
  raise ValueError("Region name required: Amygdala or sACC")
print("Running analysis on " + reg)
# define paths to data
plink_prefix_path = 'plink/MDDstudy_plink_n458'
expression_bed = 'tx_expression_'+reg+'.bed.gz'
covariates_file = 'tx_covariates_'+reg+'.txt'
prefix = 'tqtl_'+reg+'_'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covars_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

## get start time
start = time.time()

#print("\n**** Start map_cis ****")
#cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, 
#                     covariates_df=covars_df)
#tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)
#cif_df.to_csv("out_cis/"+prefix+"_cis_df_qvalues.csv", sep="\t", float_format="%.6g", compression='gzip')

print("\n**** Start map_nominal ****")

#cis.map_nominal(genotype_df, variant_df, phenotype_df,phenotype_pos_df,
#                prefix=prefix, covariates_df=covars_df,
#                output_dir= "out_nominal", verbose=False)
 
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                prefix=prefix, covariates_df=covars_df,
                output_dir= "out_nominal")

end = time.time()
 
print("runtime:", end-start)
