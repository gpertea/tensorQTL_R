#!/usr/bin/env python
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis
import time
import sys
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

freg=sys.argv[1]
pair = freg.split("_")
feature = pair[0]
region = pair[1]
print("Running tensorQTL gw nominal on " + feature+ " - " + region)

outdir='out/genomewide_nominal_062523'
tag = feature + "_" + region + '_gwnom'    # output file prefix

#expres, covar = get_input_paths(feature, region)
expres = 'tqtl_input/' + freg + '_expression.bed.gz'
covar  = 'tqtl_input/' + freg+ '_covariates.txt'
#plink = "/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"
plink = "/dcs04/lieber/lcolladotor/chessBrain_LIBD4085/mdd_old/eqtl_tf/tqtl_input/plink/MDDstudy_plink_n458"

#genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = load_data(plink, expres, covar, add_chr = True, fix_geno_names = True)
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expres)
print("Phenotype dimensions:", end='')
print(phenotype_df.shape)
covariates_df = pd.read_csv(covar, sep='\t', index_col=0).T
# PLINK reader for genotypes

#print("Reading Plink files: " + plink_prefix_path)
pr = genotypeio.PlinkReader(plink)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

print("Genotype dimensions:", end='')
print(genotype_df.shape)

### Fix genoSample names -- not needed for this use case - we have BrNums already
##print("Using fam to assign genoSample names")
##fam = pd.read_table(plink + '.fam', delimiter = " ", names = ["V" + str(i) for i in range(6)])
##genotype_df.columns = [str(v1) + "_" + str(v0) for (v0,v1) in zip(fam.V1, fam.V0)]

variant_df.chrom = [s.split(':')[0] for s in list(variant_df.index)]
# Filter expression to chrom in snp data
variant_chrom = set(variant_df.chrom)
express_chrom = set(phenotype_pos_df.chr)
shared_chrs = variant_chrom.intersection(express_chrom)
print("Chr to use: " + ' '.join(shared_chrs))
if express_chrom - variant_chrom:
    chrom_filter = phenotype_pos_df.chr.isin(variant_chrom)
    phenotype_df = phenotype_df[chrom_filter]
    phenotype_pos_df = phenotype_pos_df[chrom_filter]

print("Phenotype dimensions:", end='')
print(phenotype_df.shape)

print('\n**** STARTING tensorQTL ****')
print("Saving output to: " + outdir + '/' + tag)

cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix = tag, covariates_df= covariates_df,
                maf_threshold=0.05, window=500000, run_eigenmt=True, output_dir= outdir, write_top=False, verbose=False)

print(" DONE !")

