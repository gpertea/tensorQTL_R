#!/bin/bash

## make sure the collapsed annotation gtf was produced
# ./collapse_annotation.py --transcript_blacklist gencode24-25_unannotated_readthrough_blacklist.txt gencode25.main.gtf gencode25_collapsed.main.gtf

## also used 00_rse_gene2gct.R to generate the gct files and sample2subj.tab file
## also prepared the _genotype_pcs.txt and _additional_covariates.txt files

## also make sure the genotype file is using the proper chr prefix:
## in genotypes folder:
## plink2 --pfile ../genotypes/claustrum_maf05 --output-chr chrM --make-pgen --out ../genotypes/claustrum_maf05_chr
# cut -f1 ../genotypes/claustrum_maf05_chr.pvar | uniq > chrs.lst
tpm_gct=claustrum_gene_tpm.gct.gz
counts_gct=claustrum_gene_counts.gct.gz
annotation_gtf=gencode25_collapsed.main.gtf
sample2ind=claustrum_sample2subj.tab
## vcf_chr_list=chrs.lst ## not used, default for --chrs is to use chromosomes chr1-chr22,chrX
prefix=clau_noERCC
python eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} ${annotation_gtf} \
    ${sample2ind} ${prefix} \
    --tpm_threshold 0.1 \
    --count_threshold 6 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm

#The number of PEER factors was selected as function of sample size (N):
#- 15 factors for N < 150
#- 30 factors for 150 ≤ N < 250
#- 45 factors for 250 ≤ N < 350
#- 60 factors for N ≥ 350
num_peer=15
Rscript run_PEER.R ${prefix}.expression.bed.gz ${prefix} ${num_peer}


##This next step generates a combined covariates file, containing genotype PCs, PEER factors,
##  and additional explicit covariates (e.g., genotyping platform).

## The covariate files should have one covariate per row, with an identifier in the first column,
## and a header line with sample identifiers. This step will generate the file `${prefix}.combined_covariates.txt`
# fgrep -v ERCCsumLogErr clau_ERCC.addtl_covars.txt > clau_noERCC.addtl_covars.txt
add_covariates=${prefix}.addtl_covars.txt
genotype_pcs=claustrum_genotype_pcs.txt
python combine_covariates.py ${prefix}.PEER_covariates.txt ${prefix} \
    --genotype_pcs ${genotype_pcs} \
    --add_covariates ${add_covariates}

## run tensortqtl cis-QTL mapping with permutations

plink_prefix_path="../genotypes/claustrum_maf05_chr"

python3 -m tensorqtl ${plink_prefix_path} ${prefix}.expression.bed.gz ${prefix} \
    --covariates ${prefix}.combined_covariates.txt \
    --mode cis
