library(SummarizedExperiment)
library(sessioninfo)
library(tidyverse)
#library(VariantAnnotation)
#library(jaffelab)
library(recount)
library(data.table)
library(qs)
library(sva)

rse2bed <- function(rse, feat) {
  rr_df <- as.data.frame(rowRanges(rse))
  counts <- NULL
  if (feat=='gene') { #TODO: adjust the feature naming by feature type
     rownames(rr_df) <- paste0(rowRanges(rse)$Symbol, '|', rownames(rr_df))
     counts <- SummarizedExperiment::assays(rse)$rpkm
  } else if (feat=='tx') {
    rownames(rr_df) <- paste0(rowRanges(rse)$gene_name, '|', rownames(rr_df))
    counts <- SummarizedExperiment::assays(rse)$tpm
  }
  if (is.null(counts)) {
      counts <- SummarizedExperiment::assays(rse)$counts # must be there
      stopifnot(!is.null(counts))
      if (is.null(rowData(rse)$Length)) {
        message('Length not found in rowData, setting to 100')
        rowData(rse)$Length <- 100 ##junctions
      }
      counts <- recount::getRPKM(rse, 'Length')
  }
  rownames(counts) <- rownames(rr_df)

  rr_df$tss_start=ifelse(rr_df$strand=='-', rr_df$end, rr_df$start)
  rr_df$start <- rr_df$tss_start
  rr_df <- rr_df %>%  tibble::rownames_to_column("ID") %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::mutate(end = start + 1) %>%
    dplyr::select(`#Chr` = seqnames, start, end, ID)

  counts <- log2(counts+1)
  colnames(counts) <- rse$BrNum
  counts <- as.data.frame(counts) %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::left_join(rr_df, ., by = "ID")
  return(counts)
}

covar_format <- function(data, rn) {
  data <- as.data.frame(data)
  rownames(data) <- rn
  data <- t(data)
  data <- as.data.frame(data) %>% rownames_to_column("id")
  return(data)
}

## load RSEs
#rses <- qread('RSEs_uniqBr_by_region.qs')
## load RSEs with age 13+ samples
rses <- qread("RSEs_18plus_uniqBr_by_region.qs")
## restrict to these 3 regions for now
rses <- rses[c('DLPFC', 'HIPPO', 'Amygdala', 'Caudate')]
trse <- qread('rse_tx_wadhd_br598.qs')
#trse <- qread('rse_tx_wadhd_18plus.qs')
trses <- list()
regions=names(rses)
brnums <- c()
for (reg in regions) {
   rse <- rses[[reg]]
   brnums <- c(brnums, unique(rse$BrNum))
   fwrite(as.data.frame(rse$BrNum), file=paste0(reg, '_brnums.lst'), col.names=F)
   fwrite(data.frame(brnum=rse$BrNum, adhd=as.integer(rse$ADHD)), file=paste0(reg, '_brnums_ADHD.tab'), sep='\t', col.names=T)
   #trses[[reg]] <- trse[, colnames(rse)]
}

odir='eqtl_inputs'
if (!dir.exists(odir)) dir.create(odir)
snpPCs <- fread('genotypes/adhd_gt_n458_over18_MDS.snpPCs.csv')
snpPCs$SAMPLE_ID <- sub('^Br(\\d\\d\\d)$', 'Br0\\1', snpPCs$SAMPLE_ID)

if (length(setdiff(brnums, snpPCs$SAMPLE_ID))>0) {
  stop("Error: BrNums missing in snpPCs: ", paste(setdiff(brnums, snpPCs$SAMPLE_ID), collapse=', '))
}
setkey(snpPCs, 'SAMPLE_ID')
CUTOFF=0.15 ## for both rpkm and tpm
ffPCs <- list() ## feature PCs
modmats <- list()
modelstr='~Sex + Age  + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5'
for (reg in regions) { #filter each rse by filter cutoff
  ## remove chrX,Y,M
  rses[[reg]] <- rses[[reg]][!grepl('chr[XYM]', seqnames(rowRanges(rses[[reg]]))), ]
  trses[[reg]] <- trses[[reg]][!grepl('chr[XYM]', seqnames(rowRanges(trses[[reg]]))), ]
  ## apply expression cutoff filter
  rpkm <- assays(rses[[reg]])$rpkm
  stopifnot(!is.null(rpkm))
  tpm <- assays(trses[[reg]])$tpm
  stopifnot(!is.null(tpm))
  rses[[reg]] <- rses[[reg]][rowMeans(rpkm)>CUTOFF, ]
  trses[[reg]] <- trses[[reg]][rowMeans(tpm)>CUTOFF, ]
  rm(rpkm)
  rm(tpm)
  stopifnot(identical(colnames(rses[[reg]]), colnames(trses[[reg]])))
  ## calculate/load feature PCs
  rpd <- as.data.frame(colData(rses[[reg]]))
  ## add snpPCs
  stopifnot(identical(rses[[reg]]$BrNum, snpPCs[rpd$BrNum, ]$SAMPLE_ID))
  rpd <- cbind(rpd, snpPCs[rpd$BrNum, -1])
  mod <- model.matrix(as.formula(modelstr), data = rpd)
  modmats[[reg]] <- mod
  #message(" dim = ", dim(rpd)[1],  ' x ', dim(rpd)[2])
  frses=list(gene=rses[[reg]], tx=trses[[reg]])
  ffPCs[[reg]] <- list()
  # (gene=gffPCs[[reg]], tx=tffPCs[[reg]])
  for (feature in names(frses)) {
    fn=paste0(odir, '/', reg,'.',feature,'.featurePCs.qs')
    if (file.exists(fn)) {
      ffPCs[[reg]][[feature]] <- qread(fn)
    } else {
      stop("Error: file not found: ", fn)
      cat("Calculate", feature, "feature PCs in", reg, "..\n")
      lmx <- assays(frses[[feature]])$rpkm
      if (is.null(lmx)) {
        lmx <- assays(frses[[feature]])$tpm
      }
      lmx <- log2(lmx + 1)
      pca <- prcomp(t(lmx))
      k <- 0
      if (dim(lmx)[1]>60000) {
        k <- num.sv(lmx, mod, vfilter = 60000)
      } else {
        k <- num.sv(lmx, mod)
      }
      ffPCs[[reg]][[feature]] <- pca$x[, 1:k]
      qsave(ffPCs[[reg]][[feature]], file=fn)
      cat("   saved as",fn, "\n")
    }
  }
}

## map ffPCs rownames to BrNum
for (reg in regions) {
 frses=list(gene=rses[[reg]], tx=trses[[reg]])
 for (feature in names(frses)) {
   pc <- ffPCs[[reg]][[feature]]
   stopifnot(identical(rownames(pc), colnames(frses[[feature]])))
   #rownames(pc) <- rses[[reg]]$BrNum -- no, covar_format will do this

   rse <- frses[[feature]]
   rpd <- as.data.frame(colData(rse))
   cat(" exporting data for", reg,"feature",feature," ..\n")
   model <- modmats[[reg]][, -1]
   ## model (sample) factors and covariates
   stopifnot(identical(rownames(model), rownames(rpd)))
   cpd <- covar_format(model, rpd$BrNum)
   ## PC data, SVs (feature) covariates
   stopifnot(identical(rownames(pc), rownames(rpd)))
   cpc <- covar_format(pc, rpd$BrNum)
   ## bind model (sample) covars with feature covars (SVs) and save
   stopifnot(identical(colnames(cpd), colnames(cpc)))
   covars <- rbind(cpd, cpc)
   fwrite(covars, file = paste0(odir, '/', reg,'.', feature, '.covars.txt'),
          sep = "\t", quote = FALSE, row.names = FALSE)

   bed <- rse2bed(rse, feature)
   fwrite(bed, file = paste0(odir, '/', reg,'.', feature,'.expr.bed.gz'),
          sep = "\t", quote = FALSE, row.names = FALSE)
   rm(bed)
   cat("  ..exported.\n")
 }
}
