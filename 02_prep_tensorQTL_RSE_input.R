#!/bin/env Rscript
library(SummarizedExperiment)
library(sessioninfo)
library(tidyverse)
library(recount)
library(data.table)
library(qs)
library(sva)

###### CHANGE HERE - file prefixes and paths:
## only provide the gene RSE path, the others will be found automatically
#fgrse <- './rdata/rse_gene_n114.rda'
fgrse <- 'mdd_exprs_cutoff/rse_gene_Amygdala_n540.rda'
## path to SNP PCs file prepared with calc_snpPCs.R
#fsnp_pcs <- 'genotypes/gt_bsp12_MDS_snpPCs.csv'
#fsnp_pcs <- 'genotypes/gt_bsp12_rIDs_MDS.snpPCs.tab'
fsnp_pcs <- NULL ## NULL means they should be embedded already in colData(rse)
## note: tensorQTL script should use the matching genotypes/gt_bsp12_rIDs_MDS.bed as input!

#### model with known covariates to use in the eQTL analysis ###
modelstr='~Sex + Age  + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5'

## note that SVA will be used to calculate feature PCs (SVs)
## that are not accounted for in this base model

dsname <- 'mdd_amyg' ## dataset name prefix to use for all output files
## NOTE: change this if you change the model or any input data!
##       this is the same with the dsname variable in the python script (03_run_tensorqtl)

## file with genotype ID to RNAseq SAMPLE ID mapping, if needed
## can be set to NULL if the genoID is the same with RNA SAMPLE_ID
fgeno2rna <- NULL
#fgeno2rna <- 'genoID_to_sampleID.csv' ## 2 column file with genotype IDs and RNA sample IDs
## must be without header, 1st column: genotype IDs, 2nd column: RNA SAMPLE_IDs as in the RSEs
## if used, the same file must be given to the python script that runs tensorQTL
## NOTE: colnames of the RSEs must match the RNA SAMPLE_IDs !
#----------------------

geno2rna <- NULL

if (!file.exists(fgrse)) stop(paste("Cannot find gene RSE file", fgrse))
if (!is.null(fsnp_pcs) && !file.exists(fsnp_pcs)) stop(paste("Cannot find gene SNP PCs file", f_snp_pcs))

if (!is.null(fgeno2rna) && nchar(fgeno2rna)>1) {
  if (!file.exists(fgeno2rna)) stop(paste("Cannot find genotype ID to RNAseq SAMMPLE ID file", fgeno2rna))
  geno2rna <- fread(fgeno2rna, header=F, data.table=F)
  colnames(geno2rna) <- c('genoID', 'SAMPLE_ID')
}

fn <- basename(file.path(fgrse))
fdir <- dirname(file.path(fgrse))
fpat <- paste0(sub('_gene', '_.+', fn, ignore.case = T), '$')
frses <- list.files(fdir, pattern=fpat, full.names=T)
names(frses) <- tolower(sub('rse_([^_\\.]+).+','\\1',basename(frses)))
pref_order <- c('gene', 'exon', 'tx', 'jx', 'jxn')
ord_names <- factor(names(frses), levels = pref_order)
frses <- frses[order(ord_names)]

rse2bed <- function(rse, feat) { ## assumes tpm and rpkm assays are already present
  rr_df <- as.data.frame(rowRanges(rse))
  counts <- NULL
  if (feat=='tx') {
    if (!is.null(rr_df$gene_name) & !grepl('\\|', rownames(rr_df)[1])) {
      rownames(rr_df) <- paste0(rr_df$gene_name, '|', rownames(rr_df))
      rownames(rse) <- rownames(rr_df)
    }
    counts <- SummarizedExperiment::assays(rse)$tpm
  } else {
    if (feat %in% c('gene', 'exon') & !is.null(rr_df$Symbol)
        & !grepl('\\|', rownames(rr_df)[1])) {
        rownames(rr_df) <- paste0(rr_df$Symbol, '|', rownames(rr_df))
        rownames(rse) <- rownames(rr_df)
    }
    counts <- SummarizedExperiment::assays(rse)$rpkm
  }

  stopifnot(!is.null(counts))

  rr_df$tss_start=ifelse(rr_df$strand=='-', rr_df$end, rr_df$start)
  rr_df$start <- rr_df$tss_start
  rr_df <- rr_df %>%  tibble::rownames_to_column("ID") %>%
  dplyr::arrange(seqnames, start) %>%
  dplyr::mutate(end = start + 1) %>%
  dplyr::select(`#Chr` = seqnames, start, end, ID)
  counts <- log2(counts+1)
  colnames(counts) <- colnames(rse) # RNAseq SAMPLE_ID (was rse$BrNum)
  counts <- as.data.frame(counts) %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(rr_df, ., by = "ID")
  return(counts)
}

covar_format <- function(data) {
  data <- as.data.frame(data)
  #rownames(data) <- rn
  data <- t(data)
  data <- as.data.frame(data) %>% rownames_to_column("id")
  return(data)
}


odir='eqtl_inputs' ## output directory for expression and covariate files
if (!dir.exists(odir)) dir.create(odir)
if (!is.null(fsnp_pcs)) {
  snpPCs <- fread(fsnp_pcs, data.table=F)
  if (!is.null(geno2rna)) {
     genoIDs <- setdiff(geno2rna$genoID, snpPCs$SAMPLE_ID)
     if (length(genoIDs)>0) {
       stop("Error: genotype IDs lacking snpPCs: ", paste(genoIDs, collapse=', '))
     }
  }
  rownames(snpPCs) <- snpPCs$SAMPLE_ID ## this is the genotype ID
}
## use a basic expression cutoff
#EXPR_CUTOFF=0.2
EXPR_CUTOFF=0 ## already applied for the mdd set
pd <- NULL ## phenodata used for all RSEs
model <- NULL

## optional: reduce frses to only a subset of features
#frses <- frses[c('gene', 'tx')]
for (i in 1:length(frses)) {
  rse <- get(load(file = frses[[i]], verbose=F)) # load the RSE
  rm(list=ls(pattern='^rse_'))
  feature <- names(frses)[i]
  ## subset RSE by geno2rna$SAMPLE_ID, which should match colnames(pd)
  if (!is.null(fsnp_pcs)) {
    nc <- ncol(rse)
    if (is.null(geno2rna)) {
      rse <- rse[ ,colnames(rse) %in% snpPCs$SAMPLE_ID]
    } else {
      rse <- rse[ ,colnames(rse) %in% geno2rna$SAMPLE_ID]
    }
    if (ncol(rse)<nc/2) {
      stop("Too many RSE samples do not match the SAMPLE ID with genotype data!")
    }
  }
  if (is.null(pd)) {
    pd <- as.data.frame(colData(rse))
    if (!is.null(fsnp_pcs)) {
      if (is.null(geno2rna)) {
        ## SAMPLE_ID is the same as genotype ID
        missIDs <- setdiff(rownames(pd), rownames(snpPCs))
        if (length(missIDs)>0) {
          stop("Error: RNA SAMPLE_IDs lacking snpPCs: ", paste(missIDs, collapse=', '))
        }
        snpPCs <- snpPCs[rownames(pd), ] #subset/reorder snpPCs
        stopifnot(identical(rownames(pd), rownames(snpPCs)))
        pd <- cbind(pd, snpPCs[, grep('^snpPC\\d+', colnames(snpPCs))])
      } else {
        ## reorder geno2rna to match order in colnames(pd)
        ## append genoID to pd
        rownames(geno2rna) <- geno2rna$SAMPLE_ID
        geno2rna <- geno2rna[rownames(pd), ]
        stopifnot(identical(rownames(pd), rownames(geno2rna)))
        pd <- cbind(pd, data.frame(genoID = geno2rna$genoID))
        snpPCs <- snpPCs[pd$genoID, ]
        stopifnot(identical(pd$genoID, rownames(snpPCs)))
        pd <- cbind(pd, snpPCs[, grep('^snpPC\\d+', colnames(snpPCs))])
      }
    }
    ## model (sample) factors and covariates
    model <- model.matrix(as.formula(modelstr), data = pd) [, -1]
    stopifnot(identical(rownames(model), rownames(pd)))
  } # pd and model buing build for the first RSE
  stopifnot(identical(colnames(rse), rownames(pd)))
  if (feature=='jx' & is.null(assays(rse)$rpkm) & !is.null(assays(rse)$rp10m)) {
    names(assays(rse)) <- sub('rp10m', 'rpkm', names(assays(rse)))
  }
  if (feature=='tx') {
    stopifnot(!is.null(assays(rse)$tpm))
  } else { # not tx, must have rpkm
    if (is.null(assays(rse)$rpkm)) {
      stopifnot(!is.null(assays(rse)$counts))
      if (grepl('^jx', feature) & is.null(rowData(rse)$Length)) {
          rowData(rse)$Length <- 100 ##junctions
      }
      assays(rse)$rpkm <- recount::getRPKM(rse, 'Length')
    }
  }
  ## apply expression cutoff
  if (EXPR_CUTOFF>0) {
    if (feature=='tx') {
      rse <- rse[rowMeans(assays(rse)$tpm)>EXPR_CUTOFF, ]
    } else {
      rse <- rse[rowMeans(assays(rse)$rpkm)>EXPR_CUTOFF, ]
    }
  }
  ## discard chromosome Y and M data
  rse <- rse[!grepl('chr[YM]', seqnames(rowRanges(rse))), ]

  ## calculate feature PCs (sva)
  fn=paste0(odir, '/', dsname, '.', feature,'.featurePCs.qs')
  ffPCs <- NULL
  if (file.exists(fn)) {
    cat("Loading", feature, "feature SVs from", fn, "\n")
    ffPCs <- qread(fn)
  } else {
    cat("Calculate", feature, "feature SVs..\n")
    lmx <- assays(rse)$rpkm
    if (is.null(lmx)) {
      lmx <- assays(rse)$tpm
    }
    lmx <- log2(lmx + 1)
    pca <- prcomp(t(lmx))
    k <- 0
    if (dim(lmx)[1]>60000) {
      k <- num.sv(lmx, model, vfilter = 60000)
    } else {
      k <- num.sv(lmx, model)
    }
    ffPCs <- pca$x[, 1:k]
    qsave(ffPCs, file=fn)
    cat("     saved as",fn, "\n")
  }

  stopifnot(identical(rownames(ffPCs), colnames(rse)))

  cat(" exporting data for feature",feature," ..\n")

  cpd <- covar_format(model)
  ## PC data, SVs (feature) covariates
  stopifnot(identical(rownames(ffPCs), rownames(pd)))
  cpc <- covar_format(ffPCs)
  ## bind model (sample) covars with feature covars (SVs) and save
  stopifnot(identical(colnames(cpd), colnames(cpc)))
  covars <- rbind(cpd, cpc)
  fwrite(covars, file = paste0(odir, '/', dsname, '.', feature, '.covars.txt'),
         sep = "\t", quote = FALSE, row.names = FALSE)

  bed <- rse2bed(rse, feature)
  fbed=paste0(odir, '/', dsname, '.', feature,'.expr.bed.gz')
  fwrite(bed, file = fbed, sep = "\t", quote = FALSE, row.names = FALSE)
  rm(bed)
  cat("  ..exported to", fbed, "\n")
}

