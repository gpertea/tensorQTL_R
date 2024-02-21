library(SummarizedExperiment)
library(sessioninfo)
library(tidyverse)
#library(VariantAnnotation)
#library(jaffelab)
library(recount)
library(data.table)
library(qs)
library(sva)

###### --- change here file prefixes and paths:
## only provide the gene RSE path, the others will be found automatically
fgrse <- './rdata/rse_gene_n114.rda'

## path to SNP PCs file prepared with calc_snpPCs.R
fsnp_pcs <- 'genotypes/gt_bsp12_MDS.snpPCs.csv'

## model with known covariates to use in the eQTL analysis
modelstr='~Sex + Age  + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5'

fgeno2rna <- 'genoID_to_sampleID.csv' ## 2 column file with genotype IDs and RNA sample IDs
## must be without header, 1st column: genotype IDs, 2nd column: RNA SAMPLE_IDs as in the RSEs
## NOTE: colnames of the RSEs must match the RNA SAMPLE_IDs !
## ### --------------------------------------------

if (!file.exists(fgrse)) stop(paste("Cannot find gene RSE file", fgrse))
if (!file.exists(fsnp_pcs)) stop(paste("Cannot find gene SNP PCs file", f_snp_pcs))
if (!file.exists(fgeno2rna)) stop(paste("Cannot find genotype ID to RNAseq SAMMPLE ID file", fgeno2rna))
geno2rna <- fread(fgeno2rna, header=F, data.table=F)
colnames(geno2rna) <- c('genoID', 'SAMPLE_ID')

rses <- list()
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


pd <- NULL ## phenodata used for all RSEs
odir='eqtl_inputs' ## output directory for expression and covariate files
if (!dir.exists(odir)) dir.create(odir)

snpPCs <- fread(fsnp_pcs, data.table=F)
genoIDs <- setdiff(geno2rna$genoID, snpPCs$SAMPLE_ID)
rownames(snpPCs) <- snpPCs$SAMPLE_ID ## this is the genotype ID
if (length(genoIDs)>0) {
  stop("Error: genotype IDs lacking snpPCs: ", paste(genoIDs, collapse=', '))
}

## use a basic expression cutoff
EXPR_CUTOFF=0.2
#modmats <- list()
for (i in 1:length(frses)) {
  rse <- get(load(file = frses[[i]], verbose=F)) # load the RSE
  feature <- names(frses)[i]
  ## subset RSE by geno2rna$SAMPLE_ID, which should match colnames(pd)
  nc <- ncol(rse)
  rse <- rse[ ,colnames(rse) %in% geno2rna$SAMPLE_ID]
  if (ncol(rse)<nc/2) {
    stop(paste("Error: too many RSE samples do not match the SAMPLE ID in mapping file ",
               fgeno2rna))
  }
  if (is.null(pd)) {
    pd <- as.data.frame(colData(rse))
    ## reorder geno2rna to match order in colnames(pd)
    ## append genoID to pd
    rownames(geno2rna) <- geno2rna$SAMPLE_ID
    geno2rna <- geno2rna[rownames(pd), ]
    stopifnot(identical(rownames(pd), rownames(geno2rna)))
    pd <- cbind(pd, data.frame(genoID = geno2rna$genoID))
    snpPCs <- snpPCs[pd$genoID, ]
    stopifnot(identical(pd$genoID, rownames(snpPCs)))
    stopifnot(identical(pd$genoID, rownames(snpPCs)))
  }
  stopifnot(identical(colnames(rse), rownames(pd)))
  if (feature!='tx' & !is.null(assays(rse)$rpkm)) {
    if (grepl('^jx', feature) & is.null(rowData(rse)$Length)) {
        rowData(rse)$Length <- 100 ##junctions
    }
    assays(rse)$rpkm <- recount::getRPKM(rse, 'Length')
  }
  ## apply expression cutoff
  if (feature=='tx') {
    rse <- rses[rowMeans(assays(rse)$tpm)>EXPR_CUTOFF, ]
  } else {
    rse <- rse[rowMeans(assays(rse)$rpkm)>EXPR_CUTOFF, ]
  }
  ## discard chromosome Y and M data
  rse <- rse[!grepl('chr[YM]', seqnames(rowRanges(rse))), ]

  ## calculate feature PCs (sva)
  fn=paste0(odir, '/', feature,'.featurePCs.qs')
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
      k <- num.sv(lmx, mod, vfilter = 60000)
    } else {
      k <- num.sv(lmx, mod)
    }
    ffPCs <- pca$x[, 1:k]
    qsave(ffPCs, file=fn)
    cat("     saved as",fn, "\n")
  }
  pc <- ffPCs
  stopifnot(identical(rownames(pc), colnames(rse)))
  cat(" exporting data for feature",feature," ..\n")
  model <- model.matrix(as.formula(modelstr), data = pd) [, -1]
  ## model (sample) factors and covariates
  stopifnot(identical(rownames(model), rownames(pd)))

  cpd <- covar_format(model)
  ## PC data, SVs (feature) covariates
  stopifnot(identical(rownames(pc), rownames(pd)))
  cpc <- covar_format(pc)
  ## bind model (sample) covars with feature covars (SVs) and save
  stopifnot(identical(colnames(cpd), colnames(cpc)))
  covars <- rbind(cpd, cpc)
  fwrite(covars, file = paste0(odir, '/', reg,'.', feature, '.covars.txt'),
         sep = "\t", quote = FALSE, row.names = FALSE)

  bed <- rse2bed(rse, feature)
  fbed=paste0(odir, '/',feature,'.expr.bed.gz')
  fwrite(bed, file = fbed, sep = "\t", quote = FALSE, row.names = FALSE)
  rm(bed)
  cat("  ..exported to", fbed, "\n")
}

