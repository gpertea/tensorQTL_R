library(SummarizedExperiment)
library(sessioninfo)
library(tidyverse)
library(VariantAnnotation)
library(jaffelab)
library(here)

#source(here("eqtl", "code", "rse_to_bed.R"))
## actually defined in utils.R:
read_adj_filter <- function(parquet_files, cutoff = 0.05, pval_name = "pval_nominal") {
  eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) %>%
    mutate(FDR = p.adjust(get(pval_name), "fdr"))
  message("n pairs: ", nrow(eqtl_out))
  # filter
  eqtl_out <- eqtl_out %>%
    filter(FDR < cutoff)
  message("n pairs FDR<", cutoff, ": ", nrow(eqtl_out))
  return(eqtl_out)
}

rse_to_bed <- function(rse) {
  rr_df <- as.data.frame(rowRanges(rse))
  rr_df$tss_start=ifelse(rr_df$strand=='-', rr_df$end, rr_df$start)
  rr_df$start <- rr_df$tss_start
  rr_df <- rr_df %>%  tibble::rownames_to_column("ID") %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::mutate(end = start + 1) %>%
    dplyr::select(`#Chr` = seqnames, start, end, ID)

  #counts <- SummarizedExperiment::assays(rse)$counts
  counts <- SummarizedExperiment::assays(rse)$counts
  if (!is.null(SummarizedExperiment::assays(rse)$tpm)) {
    counts <- SummarizedExperiment::assays(rse)$tpm
  } else {
    if (is.null(rowData(rse)$Length)) {
      rowData(rse)$Length <- 100
    }
    counts <- recount::getRPKM(rse, 'Length')
  }
  counts <- log2(counts+1)

  colnames(counts) <- rse$genoSample
  # counts = round(counts, 3) ## round to 3 ?

  counts <- as.data.frame(counts) %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::left_join(rr_df, ., by = "ID")

  return(counts)
}

## load
load(here("mdd_old/exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("mdd_old/exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
load(here("mdd_old/exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
load(here("mdd_old/exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## Split gene data
regions <- c(amyg = "Amygdala", sacc = "sACC")
features <- c("gene", "exon", "jxn", "tx")
names(features) <- features

rse_gene_split <- map(regions, ~ rse_gene[, rse_gene$BrainRegion == .x])
samples_split <- map(rse_gene_split, colnames)

#### Covariate Data ####
pcs <- map(features, function(f_name) map(regions, ~ get(load(here("mdd_old/eqtl", "data", "featuresPCs", paste0(f_name, "PCs", .x, ".rda"))))))
corner(pcs$gene$amyg)

covar_format <- function(data, rn) {
    data <- as.data.frame(data)
    rownames(data) <- rn
    data <- t(data)
    data <- as.data.frame(data) %>% rownames_to_column("id")
    return(data)
}

covars <- map2(
    pcs, features,
    function(pc, feat) {
        pmap(list(rse = rse_gene_split, region = regions, pc = pc), function(rse, region, pc) {
            message(paste(feat, region))
            ## Phenodata
            pd <- as.data.frame(colData(rse)[, c("PrimaryDx", "Sex", paste0("snpPC", 1:5))])
            pd <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)[, 2:9]
            pd <- covar_format(pd, rse$genoSample)

            ## PC data
            pc <- covar_format(pc, rse$genoSample)
            ## bind and save
            covars <- rbind(pd, pc)
            write.table(covars,
                file = here("mdd_old/eqtl", "gdata", "covariates_txt", paste0("covariates_", feat, "_", region, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE
            )
            return(covars)
        })
    }
)
corner(covars$gene$amyg)

#### Expression Data ####
expression_fn <- map(features, function(feat) map(regions,
#                                ~ here("mdd_old/eqtl", "gdata", "expression_bed", paste0(feat, "_", .x, ".bed"))))
                        ~ here("mdd_old/eqtl", "gdata", "expression_bed", paste0(feat, "_", .x, ".bed.gz"))))

expression_bed <- map2(list(rse_gene, rse_exon, rse_jxn, rse_tx), features, function(rse, feat) {
    rse_split <- map(regions, ~ rse[, rse$BrainRegion == .x])
    message("generating expression BED for ",feat, "..")
    expr_bed <- map(rse_split, rse_to_bed)
    message("    ", feat, " BED done.")
    return(expr_bed)
})

walk2(expression_bed, expression_fn, function(expr, fn) {
    #walk2(expr, fn, ~ write.table(.x, .y,
    walk2(expr, fn, ~ data.table::fwrite(.x, file=.y,
        sep = "\t",
        quote = FALSE, row.names = FALSE
    ))
})

## Create shell script to zip data
#commands <- map(unlist(expression_fn), ~ paste0("bgzip ", .x, " && tabix -p bed ", .x, ".gz"))
#if (file.exists("bed_bgzip.sh")) file.remove("bed_bgzip.sh")
#sgejobs::job_single("bed_bgzip",
#    create_shell = TRUE, queue = "bluejay", memory = "100G",
#    command = paste(commands, collapse = "\n")
#)
## Add "module load htslib"

#### VCF ####
#risk_vcf <- readVcf(here("eqtl", "data", "risk_snps", "LIBD_maf01_gwas_BPD.vcf.gz"))
#risk_vcf
#risk_vcf_split <- map(rse_gene_split, ~ risk_vcf[, .x$genoSample])
#map(risk_vcf_split, dim)

#vcf_fn <- map(regions, ~ here("eqtl", "data", "risk_snps", paste0("LIBD_maf01_gwas_BPD_", .x, ".vcf.gz")))
#walk2(risk_vcf_split, vcf_fn, ~ writeVcf(.x, .y))

## plink commands
#map(vcf_fn, ~ paste("plink --make-bed --output-chr chrM --vcf", .x, "--out", gsub(".vcf.gz", "", .x)))


## check

#map2(bed, risk_vcf_split, ~ all(colnames(.x[, 5:ncol(.x)]) == colnames(.y)))
#map2(bed, covars, ~ all(colnames(.x[, 5:ncol(.x)]) == colnames(.y[[1]][, 2:ncol(.y[[1]])])))


## prep interaction csv
walk2(rse_gene_split, regions, function(rse, region) {
    cell_fractions <- colData(rse)[, c("Astro", "Endo", "Macro", "Micro", "Mural", "Oligo", "OPC", "Tcell", "Excit", "Inhib")]
    cell_fractions <- as.data.frame(cell_fractions)
    rownames(cell_fractions) <- rse$genoSample
    write.csv(cell_fractions, file = here("mdd_old/eqtl", "gdata", "interaction", paste0("cell_fraction_", region, ".csv")))
})

## create shell commands ##
#sgejobs::job_single("tensorqtl_risk_snps",
#    create_shell = TRUE, queue = "bluejay", memory = "50G",
#    command = "python tensorqtl_risk_snps.py"
#)
