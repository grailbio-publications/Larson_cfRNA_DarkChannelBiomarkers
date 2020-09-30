#' Non-parametric Wilcoxon significance test between groups
#'
#' @param rpm_df dataframe with columns Subtype, gene_id,
#'  and expression_new
#' @param gene character gene name
#' @param group1 first group in dataframe to test
#' @param group2 second group in dataframe to test
#'
#' @return a numeric p-value
#' @export
#'
wilcoxon_test <- function(rpm_df, gene, group1, group2) {
  gene_dat <- rpm_df %>%
    dplyr::mutate(Subtype = as.character(Subtype)) %>%
    dplyr::filter(gene_id == gene,
                  Subtype == group1 |
                    Subtype == group2)
  return(stats::wilcox.test(expression_new ~ Subtype, data=gene_dat)$p.value)
}
