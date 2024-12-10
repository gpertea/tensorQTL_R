args <- commandArgs(trailingOnly = TRUE)
usage_msg  <- "Usage: Rscript 00_rse_gene2gct.R <gene_rse_file_path> <dataset_name>"
if (length(args) < 2) {
  stop(usage_msg)
}

frse <- args[1]
dsname <- args[2]

if (!file.exists(frse)) {
  stop("Error: The specified gene RSE file does not exist.")
}

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
})

# gene_file <- rse_files[1]
genv <- new.env()
load(file.path("data", frse), envir = genv)
if (exists("rse_gene", envir = genv)) {
    rse <- genv$rse_gene
} else {
    stop("No rse_gene object found in the gene file")
}


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
