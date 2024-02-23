#!/bin/env Rscript
## script to get MDS and snpPCs for a set of genotypes in a VCF.gz input data
## on JHPCE, the plink module must be loaded (plink must be available in PATH)!
library(data.table)
## CHANGE HERE the path to the genotype file to calculate MDS for:
#gtinput='genotypes/gt_bsp12_n114.vcf.gz' ## if you have a vcf
gtinput='genotypes/gt_bsp12_n114_rIDs.bed' ## if you have a plink bed file instead of vcf

## if you have a plink file instead of vcf, that you want to QC and use for MDS calculation
## set gtinput to the .bed file directly (NOTE: include the .bed extension here!) but
## set bedout to a different prefix! (to prevent overwriting)

bedout="genotypes/gt_bsp12_rIDs_MDS"
## write the final/filtered bed files with this bedout prefix

## QC filter
qcflt='--geno 0.01 --maf 0.001 --hwe 0.0001'

## this should probably not be applied to an already extracted VCF that has all the samples of interest
## applying the maf QC filter for very small sets may remove too many variants
## e.g. for a file with 68 genotypes, using qcflt like this can eliminate up to a half of the variants!
plink2=Sys.which('plink2')
if (nchar(plink2)==0) stop("plink2 not found in PATH!")

plink=Sys.which('plink')
if (nchar(plink)==0) stop("plink not found in PATH!")

## gtinput file must exists
if (!file.exists(gtinput)) stop(paste("Cannot find input genotype file:", gtinput))

is.vcf=TRUE
in.fmt=paste('--vcf', gtinput)
if (grepl('\\.bed$', gtinput)) {
  is.vcf=FALSE
  gtinput=sub('\\.bed$', '', gtinput)
  in.fmt=paste('--bfile', gtinput)
}

## convert vcf to plink BED (or just QC filter the input BED if given directly)
## NOTE: remove qcflt if this is a very small sub-sample already
cmd=paste(plink2,"--threads 8 --make-bed --not-chr X --output-chr chrM",
          in.fmt, qcflt, '--out',   bedout)
message(" Running command:\n",cmd)
  ## add '--keep samples_to_extract.txt' if only a subset of genotyping samples are needed
r=system(cmd)
if (r!=0) stop("Error in plink2 command!")

indfile=paste0(bedout, '_indep')

### independent and cluster
## --indep <window size>['kb'] <step size (variant ct)> <VIF threshold>
## produce a pruned subset of markers that are in approximate linkage equilibrium with each other
## --indep requires three parameters:
##          a window size in variant count or kilobase (if the 'kb' modifier is present) units,
##          a variant count to shift the window at the end of each step,
##          a variance inflation factor (VIF) threshold
cmd = paste(plink, "--bfile", bedout, "--indep 100 10 1.25 --out", bedout)
system(cmd)

## MDS components
## use "--mds-plot 15" below if you want 15 SNPs
cmd = paste0(plink, " --bfile ", bedout,
             " --cluster --mds-plot 10 --extract ", bedout, ".prune.in --out ", bedout)
system(cmd)

# ## A transpose
cmd = paste(plink, "--bfile", bedout, "--recode A-transpose --out", bedout)
system(cmd)

#### read in MDS
mds = read.table(paste0(bedout, ".mds"), header=TRUE,as.is=TRUE)
rmds = mds[,-(1:3)] #remove FID, IID and SOL columns 1-3
colnames(rmds) = paste0("snpPC",1:ncol(rmds))

## WARNING: the assumption here is that FID does not matter (it is not used)
rmds$SAMPLE_ID=mds$IID

setcolorder(rmds, 'SAMPLE_ID')

## write snpPCs file as csv
fwrite(rmds, file=paste0(bedout, ".snpPCs.tab"), sep='\t')
