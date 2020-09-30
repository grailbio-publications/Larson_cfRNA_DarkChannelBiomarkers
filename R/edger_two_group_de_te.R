#' Find differentiated expressed genes in two classes of RNA-seq data
#'
#' Takes in a vector of control / housekeeping genes to use for expression
#' analysis of a targeted RNA panel. edgeR normalization assumes that most of
#' the genes are not differentially expressed, which might not be true in a
#' targeted panel (https://support.bioconductor.org/p/87121/).
#'
#' @param genetable a data_frame with columns gene_id and columns corresponding
#' to the numeric value of RNA counts per gene in a sample.
#' @param labels a vector of length ncol(genetable)-1, corresponding to the sample groups
#' @param controls vector of genes that you expect to have consistent expression between groups
#' @param method character, normalization method used by edgeR
#'
#' @return data frame with differential p-values for each gene
#' @export
#'
edger_two_group_de_te <- function (genetable, labels, controls,
                                   method = 'TMMwsp') {
  counts <- genetable_to_matrix(genetable, transpose = F)
  label_table <- data.frame(sample_id = colnames(counts), class = factor(labels))
  # design matrix
  design <- stats::model.matrix(~0 + class, data = label_table)
  d <- edgeR::DGEList(count = counts, group = labels)
  # estimate normalization factors
  dcon <- d[controls, ]
  dcon <- edgeR::calcNormFactors(dcon, method = method)
  d$samples$norm.factors <- dcon$samples$norm.factors
  # estimate overall dispersion for dataset (generates common.dispresion and AveLogCPM)
  d <- edgeR::estimateGLMCommonDisp(d, design, verbose = TRUE)
  # estimate gene-wise dispersion estimates (generates trend.methoc and trend.dispersion)
  d <- edgeR::estimateGLMTrendedDisp(d, design)
  # estimate gene-wise dispersion estimates (generates span, prior.df, tagwise.dispersion)
  d <- edgeR::estimateGLMTagwiseDisp(d, design)
  # fit glm
  fit <- edgeR::glmFit(d, design)
  # likelihood ratio tests for coefficients
  lrt.term <- edgeR::glmLRT(fit, contrast = c(1, -1))
  # get table with pval, FDR, LR
  top.term <- edgeR::topTags(lrt.term, n = nrow(counts))
  gene_list.top.term <- as.data.frame(top.term)
  RNAseq_pvalue_biomarker <- gene_list.top.term[match(rownames(counts),
                                                      rownames(gene_list.top.term)), ]
  RNAseq_pvalue_biomarker <- RNAseq_pvalue_biomarker %>%
    dplyr::mutate(gene_id = row.names(.)) %>%
    dplyr::select(gene_id, dplyr::everything())
  return(RNAseq_pvalue_biomarker)
}
