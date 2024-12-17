args <- commandArgs(trailingOnly = TRUE)
usage_msg  <- "Usage: Rscript 00_rse_gene2gct.R <gene_rse_file_path> <dataset_name>"
if (length(args) < 2) {
  stop(usage_msg)
}

frse <- args[1]
dsname <- args[2]

## inline test
frse <- "../data/rse_gene_cla_new_n60.qs"
dsname <- "claustrum"

if (!file.exists(frse)) {
  stop("Error: The specified gene RSE file does not exist.")
}

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(qs)
})

# gene_file <- rse_files[1]
if (grepl('.Rdata$', frse)) {
  genv <- new.env()
  load(file.path("data", frse), envir = genv)
  if (exists("rse_gene", envir = genv)) {
      rse <- genv$rse_gene
  } else {
      stop("No rse_gene object found in the gene file")
  }
} else { # must be .qs
  rse <- qread(frse)
}

pd <- as.data.frame(colData(rse))
dxset  <- unique(pd$Dx)
if ('Control' %in% dxset) {
  pd$Dx <- factor(as.character(pd$Dx), levels = c('Control', setdiff(dxset, 'Control')))
}

modelStr='~Dx + Sex + Age + ERCCsumLogErr'
designMat <- model.matrix(as.formula(modelStr), data = pd)
# Transpose the design matrix and convert to data frame
transDesignMat <- t(designMat)
transDesignDf <- as.data.frame(transDesignMat)
stopifnot(identical(rse$SAMPLE_ID, colnames(transDesignDf)))
colnames(transDesignDf) <- rse$BrNum ## have to switch to subject ID
# Add the 'ID' column
transDesignDf <- cbind(ID = rownames(transDesignDf), transDesignDf)

# Write to file in the expected format
fwrite(transDesignDf, paste0(dsname, "_additional_covariates.txt"), row.names = FALSE, sep = "\t")


get_TPM <- function(counts, lengths) {
  # Convert lengths from base pairs to kilobases
  lengths_kb <- lengths / 1000
  # Calculate reads per kilobase (RPK)
  rpk <- counts / lengths_kb
  # Calculate scaling factors (per million RPK)
  scaling_factors <- colSums(rpk)
  # Compute TPM
  tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6
  return(tpm)
}

# Function to write GCT files
write_gct <- function(filename, counts, gene_names) {
  # Open a gzipped connection for writing text
  con <- gzfile(filename, "wt")
  # Write the two header lines required by GCT format
  writeLines("#1.2", con)
  writeLines(paste(nrow(counts), ncol(counts) + 1, sep = "\t"), con) # +1 for Description column
  # Prepare the data.table with Name, Description, and counts
  df <- data.frame(Name = rownames(counts),
                   Description = gene_names,
                   counts,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)
  # Write data
  write.table(
    df,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  close(con)
}

sample2subj <- as.data.frame(colData(rse)[, c("SAMPLE_ID", "BrNum")])
rownames(sample2subj) <- colData(rse)$BrNum
colnames(sample2subj) <- c("sample_id", "participant_id")

fwrite(sample2subj, paste0(dsname, "_sample2subj.tab"), row.names = FALSE, col.names = FALSE, sep = "\t")

gtPcs <- fread(paste0('../genotypes/', dsname, '_maf05_snpPCs.tab'), data.table = FALSE)
selPcs <- gtPcs[, paste0('snpPC',seq(1,3))] ## only the first 3 PCs
rownames(selPcs) <- gtPcs$BrNum
## we have to make sure the order matches the one in the rse
selPcs <- selPcs[rownames(sample2subj), ]
colnames(selPcs) <- c("snpPC1", "snpPC2", "snpPC3")
transPcs <- t(selPcs)
rownames(transPcs) <- colnames(selPcs)
# Convert to data frame and add the sample IDs as column names
transPcs <- as.data.frame(transPcs)
colnames(transPcs) <- gtPcs$BrNum
# Add the 'ID' column
transPcs <- cbind(ID = rownames(transPcs), transPcs)
fwrite(transPcs, paste0(dsname, "_genotype_pcs.txt"), row.names = FALSE, sep = "\t")


fgct_tpm <- paste0(dsname, "_gene_tpm.gct.gz")
fgct_cnt <- paste0(dsname, "_gene_counts.gct.gz")

gene_names <- rowData(rse)$Symbol
write_gct(fgct_cnt, assays(rse)$counts, gene_names)

tpm <- get_TPM(assays(rse)$counts, rowData(rse)$Length)
write_gct(fgct_tpm, tpm, gene_names)

