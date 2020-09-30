#' Convert the gene_id to gene_symbols
#'
#' @param data a data frame for the expression matrix. The first column is the
#' feature name and each column corresponds to a sample.
#'
#' @return the data whose first column is converted to gene symbol.
#' @export
#'
convert_gene_id <- function(data) {
  # Read the gene_id/ gene_symbol conversion data frame.
  gene2symbol <- readRDS(system.file("extdata", "gene2symbol_v19.rds",
                                     package = "cellfreetranscriptome"))

  data <- data %>% dplyr::left_join(gene2symbol, by="gene_id")
  data$gene_id <- data$gene_name
  data$gene_id[is.na(data$gene_id)] <- 'unknown'
  data$gene_name <- NULL
  data$gene_id <- make.unique(as.character(data$gene_id), sep = ".")
  return(data)
}
