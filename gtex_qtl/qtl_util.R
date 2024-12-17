suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
})

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


library(Rcpp)

cppFunction("
    NumericMatrix calc_rpkm_rcpp(NumericMatrix counts, NumericVector lengths_kb,
                  NumericVector lib_million) {
      int nrow = counts.nrow();
      int ncol = counts.ncol();
      NumericMatrix rpkm(nrow, ncol);
      for(int i = 0; i < nrow; ++i){
        for(int j = 0; j < ncol; ++j){
          rpkm(i,j) = counts(i,j) / lengths_kb[i] / lib_million[j];
        }
      }
      return rpkm;
    }
  ")


get_RPKM <- function(counts, lengths, lib_sizes = colSums(counts)) {
  counts <- as.matrix(counts)
  lengths_kb <- as.numeric(lengths) / 1000
  lib_million <- as.numeric(lib_sizes) / 1e6
  rpkm <- calc_rpkm_rcpp(counts, lengths_kb, lib_million)
  dimnames(rpkm) <- dimnames(counts)
  return(rpkm)
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
