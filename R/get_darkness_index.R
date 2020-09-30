#' Get gene darkness index from pilot healthy
#'
#' @param spliced_count_file a data frame for the expression matrix. The first column is the
#' feature name and each column corresponds to a sample.
#'
#' @return data frame containing the median and standard deviation expression values for each gene
#' @export
#'
get_darkness_index <- function(spliced_count_file) {
  count_cfrna_df <- readr::read_tsv(spliced_count_file)
  noncancer_counts <- count_cfrna_df %>%
    dplyr::select(gene_id, dplyr::starts_with('NC'))
  # sample RPMs
  rpm_nc <- get_rpm(noncancer_counts)
  # darkness index
  median_healthy <- apply(rpm_nc %>%
                            dplyr::select(-gene_id),
                          1,
                          stats::median)
  sd_healthy <- apply(rpm_nc %>%
                        dplyr::select(-gene_id),
                      1,
                      stats::sd)
  darkness_index <- data.frame(gene = rpm_nc$gene_id,
                               median = median_healthy,
                               sd = sd_healthy,
                               stringsAsFactors = FALSE)
  return(darkness_index)
}
