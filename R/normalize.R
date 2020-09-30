#' RNA-seq reads normalization using size factor
#'
#' The definition of size factor is in this paper:
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106}
#'
#' @param count_raw : The numeric matrix of raw reads count.
#' @param sample_anno : The metadata of the samples.
#'
#' @return The normalized reads count matrix.
#' @export
#'
count_normalization <- function(count_raw, sample_anno) {
  count_htseq2 <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_raw,
    colData = sample_anno,
    design = ~ 1)
  count_htseq2 <- DESeq2::estimateSizeFactors(count_htseq2)
  sample_anno$size_factor <- count_htseq2@colData$sizeFactor
  count_norm <- round(t(t(count_raw) / sample_anno$size_factor), digits = 0)
  return(count_norm)
}

#' Convert a gene expression data table into a caret-friendly matrix
#'
#' Gene expression data tables follow tidy conventions where a separate gene_id column is used in
#' place of row names. Caret expects explanatory variables (in this case genes) to be columns.
#' This function transposes the rna expression table and adds the gene_ids as column names.
#'
#' @param expression A table with a column "gene_id" that contains expression estimates.
#' @param transpose Set to \code{TRUE} to transpose rows and columns.
#'
#' @return A matrix where genes are columns and samples are rows.
#' @export
#'
genetable_to_matrix <- function(expression, transpose=TRUE) {
  gene_ids <- expression$gene_id
  expression <- expression %>%
    dplyr::select(-gene_id) %>%
    data.matrix()

  if (transpose) {
    expression <- t(expression)
    colnames(expression) <- gene_ids
  } else {
    rownames(expression) <- gene_ids
  }

  expression
}

#' Convert a caret-friendly matrix into a gene expression data table
#'
#' This function transposes the rna expression table and adds the gene_ids column.
#' Gene expression data tables follow tidy conventions where a separate gene_id column is used in
#' place of row names. Caret expects explanatory variables (in this case genes) to be columns.
#'
#' @param expression_matrix A matrix where column names correspond to gene_id.
#' @param genes_as_cols In classical ML, features are columns. Sometimes we have matrices where
#' features (genes) are rows. Set this to \code{FALSE} in order to handle this case.
#'
#' @return A data table with gene_id as a column.
#' @export
#'
matrix_to_genetable <- function(expression_matrix, genes_as_cols = TRUE) {
  if (genes_as_cols) {
    gene_ids <- colnames(expression_matrix)
    expression <- tibble::as_tibble(t(expression_matrix))
  } else {
    gene_ids <- rownames(expression_matrix)
    expression <- tibble::as_tibble(expression_matrix)
  }
  expression$gene_id <- gene_ids
  expression %>%
    dplyr::select(gene_id, dplyr::everything())
}

#' Normalize the raw count data by size factor
#'
#' @param genetable A data frame for the expression matrix. The first column is the
#' feature name and each column corresponds to a sample.
#'
#' @return A data frame with rpm data.
#' @export
#'
get_rpm <- function(genetable) {
  counts <- genetable_to_matrix(genetable, transpose=F)
  counts <- scale(counts, center=FALSE, scale=colSums(counts)) * 1e6
  matrix_to_genetable(counts, genes_as_cols=F)
}

#' Normalize reference matrix for tissue deconvolution
#'
#' @param signatures a data_frame reference matrix
#'
#' @return normalized reference matrix
#' @export
#'
normalize_matrix <- function(signatures) {
  signatures_norm <- signatures
  signatures_sum <- apply(signatures, 2, sum)
  signatures_median <- stats::median(signatures_sum)
  for (n in seq_len(dim(signatures_norm)[2])) {
    signatures_norm[, n] <- round((signatures_norm[, n]/signatures_sum[n]),
                                  digits = 3)
  }
  return(signatures_norm)
}
