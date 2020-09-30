#' Summarises gene-wise RPMs across a vector of samples
#'
#' @param gene_ids vector of gene ids
#' @param rpms matrix of gene-wise rpms, number of rows are length(gene_ids)
#' @param th_rpm numeric threshold RPM limit to binarise gene as detected (1) or not detected (0)
#'
#' @return A data_frame with one row per gene id, columns containing number detected above th_rpm,
#' sum_rpms by gene, mean_rpm by gene, and sd_rpm (standard deviation) by gene
#' @export
#'
summarise_rpms <- function(gene_ids, rpms, th_rpm = 0) {
  nums <- apply(rpms, 1, function(x) sum(x > th_rpm))
  sums <- apply(rpms, 1, function(x) sum(x))
  means <- apply(rpms, 1, function(x) mean(x))
  sds <- apply(rpms, 1, function(x) stats::sd(x))
  cbind.data.frame(gene_id = gene_ids,
                   number_detected = nums,
                   sum_rpm = sums,
                   mean_rpm = means,
                   sd_rpm = sds)
}
