aggregate_sample_counts <- function(ocx_files) {
  tables <- lapply(ocx_files, function(fp) {
    sample_id = sub("PCMP_", "", basename(fp))
    sample_id = sub("_assembled.fastq.taxa.tsv", "", sample_id)
    table <- subset(read.delim(fp), select=c(tax_id, readcount))
    table[[sample_id]] <- table$readcount
    table$readcount <- NULL
    table
  })
  Reduce(tables, f = function(...) merge(..., by="tax_id", all=TRUE))
}

aggregate_taxonomy <- function(ocx_files) {
  tables <- lapply(ocx_files, function(fp) {
    subset(read.delim(fp), select=-readcount)
  })
  Reduce(tables, f = function(...) merge(..., all=TRUE))
}

## ---- Called from Snakefile ----
if (exists("snakemake")) {
  samples <- aggregate_sample_counts(snakemake@input[[1]])
  taxa <- aggregate_taxonomy(snakemake@input[[1]])
  write.table(samples, file=snakemake@output[['all_samples']], quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(taxa, file=snakemake@output[['all_taxa']], quote = FALSE, row.names = FALSE, sep="\t")
}

