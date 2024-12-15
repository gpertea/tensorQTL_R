#!!/usr/bin/env Rscript
## get expression PCs from normalized, expression-cutoff data
library(qs)
library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(recount)
#library(jaffelab)
library(sva)
library(purrr)

setDTthreads(4)

### CHANGE rse_gene path here as needed
fgrse <- './rdata/rse_gene_n114.rda'
### -------------------------

rses <- list()
fn <- basename(file.path(fgrse))
fdir <- dirname(file.path(fgrse))
fpat <- paste0(sub('_gene', '_.+', fn, ignore.case = T), '$')
frses <- list.files(fdir, pattern=fpat, full.names=T)
rses <- map(frses, ~get(load(file = .x, verbose=T)))
names(rses) <- tolower(sub('rse_([^_\\.]+).+','\\1',basename(frses)))

pd <- as.data.frame(colData(rses[['gene']]))

#### Model ####
pd$Dx <- as.factor(as.character(pd$Dx))
pd$Dx <- relevel(pd$Dx, ref='Control')
pd$Region <- as.factor(as.character(pd$Region))
pd$BrNum <- sub('^Br0', 'Br', pd$BrNum)

### Change the model here:
modelstr=paste0('~Dx + Age + Sex + ', paste0(grep('^snpPC', colnames(pd), value=T), collapse='+') )
###
### IMPORTANT: BrNum will be used as the sample ID in the output files
model <- model.matrix(as.formula(modelstr), data = pd)

#jaffelab::expression_cutoff(assays(rses[['tx']])$tpm)
#2023-05-21 18:55:36 the suggested expression cutoff is 0.2
#percent_features_cut  samples_nonzero_cut
#0.24                 0.15
CUTOFF=0.1
## assumes a 'tpm' or 'rpkm' assay is available for all features ()
regions=unique(as.character(pd$Region))
features=names(rses)
xpms <- list()
for (feat in features) {
  #message("for feature: ", feat)
  anames=names(assays(rses[[feat]]))
  an=grep('^[tr]p', anames, value=T)
  if (length(an)==0) stop(paste("Cannot find tpm or rp*m normalized assay in", feat, 'RSE'))
  xpm <- assays(rses[[feat]])[[an]]
  colData(rses[[feat]]) <- DataFrame(pd)
  rses[[feat]] <- rses[[feat]][rowMeans(xpm)>CUTOFF, ]
  xpms[[feat]] <- an  # store assay name per RSE in xpms
  message(feat, " RSE filtered by expression cutoff ", CUTOFF, " from ",
          dim(xpm)[1], " to ", dim(rses[[feat]])[1], " features.")
  ## remove chr Y and M
  rses[[feat]] <- rses[[feat]][!grepl('chr[YM]', seqnames(rowRanges(rses[[feat]]))), ]
  rm(xpm)
}

## remove mappings on chrY, chrM

## when using risk SNPS: there are no risk snps on chr21, chrX,
#rse_tx <- rse_tx[! seqnames(rowRanges(rse_tx)) %in% c('chr21', 'chrM', 'chrX', 'chrY'), ]

odir='/eqtl'
if (!dir.exists(odir)) dir.create(odir)

regpd <- split(pd, pd$Region)
ffPCs <- list()
for (reg in regions) {
  rpd <- regpd[[reg]]
  model <- model.matrix(as.formula(modelstr), data = rpd)
  #message(" dim = ", dim(rpd)[1],  ' x ', dim(rpd)[2])
  for (feat in features) {
    ## calculating featurePCs can take a while, see if we already saved them
    freg=paste0(feat,'_',reg)
    fn=paste0(odir, '/', freg,'_featurePCs.qs')
    if (file.exists(fn)) {
      ffPCs[[freg]] <- qread(fn)
    } else {
      cat("Calculate feature PCs for",feat, "in", reg, "..\n")
      lmx <- assays(rses[[feat]])[[xpms[[feat]]]][, rpd$SAMPLE_ID]
      lmx <- log2(lmx + 1)
      mod_region <- model[rpd$SAMPLE_ID, ]
      pca <- prcomp(t(lmx))
      k <- 0
      if (dim(lmx)[1]>60000) {
        k <- num.sv(lmx, mod_region, vfilter = 60000)
      } else {
        k <- num.sv(lmx, mod_region)
      }
      ffPCs[[freg]] <- pca$x[, 1:k]
      qsave(ffPCs[[freg]], file=fn)
      cat("   saved as",fn, "\n")
    }

  }
}

covar_format <- function(data, rn) {
  data <- as.data.frame(data)
  rownames(data) <- rn
  data <- as.data.frame(t(data))
  setDT(data, keep.rownames = 'id')
  #setcolorder(data, 'id')
  return(data)
}

rse2bed <- function(rse, feat) {
  gr <- rowRanges(rse)
  if (feat=='gene') {
    names(gr) <- gr$gene_id
  } else {
    names(gr) <- gr$transcript_id
  }
  gr$feat_id <- names(gr)

  rr_df <- as.data.frame(gr)
  rr_df$start <- fifelse(rr_df$strand=='-', rr_df$end, rr_df$start)
  rr_df <- rr_df %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::mutate(end = start + 1) %>%
    dplyr::select(`#Chr` = seqnames, start, end, ID)

  counts <- log2(assays(rse)[[xpms[[feat]]]]+1)
  colnames(counts) <- colData(rse)$BrNum
  rownames(counts) <- gr$feat_id

  counts <- as.data.frame(counts) %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::left_join(rr_df, ., by = "ID")
  return(counts)
}

## now save expression data (RSE to BED) per region
## using BrNums instead of SAMPLE_ID !
#brnums <- fread('mdd_data/genotypes/brnum_samples.lst', header=F)
for (reg in regions) {
  for (feat in features) {
    rpd <- regpd[[reg]]
    rse <- rses[[feat]][, rpd$SAMPLE_ID]
    cat(" exporting data for ", feat, '_', reg, " ..\n", sep='')

    ## Phenodata
    #cpd <- rpd[, c("Dx", "Sex", paste0("snpPC", 1:5))]
    model <- model.matrix(as.formula(modelstr), data = rpd) [, -1]
    cpd <- covar_format(model, rpd$BrNum)
    ## PC data
    freg=paste0(feat, '_', reg)
    #fn=paste0(odir, '/', freg,'_featurePCs.qs')
    #pc <- qread(fn)
    pc <- ffPCs[[freg]]
    cpc <- covar_format(pc, rpd$BrNum)
    ## bind and save
    covars <- rbind(cpd, cpc)
    fwrite(covars, file = paste0(odir, '/', freg,'_covariates.txt'),
                sep = "\t", quote = FALSE, row.names = FALSE)

    bed <- rse2bed(rse, feat)
    fwrite(bed, file = paste0(odir, '/', freg, '_expression.bed.gz'),
           sep = "\t", quote = FALSE, row.names = FALSE)
    rm(bed)
    cat("  ..exported.\n")
    }
}

