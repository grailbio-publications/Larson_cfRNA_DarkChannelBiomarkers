#' Translate p-values from significance testing to astericks
#'
#' @param p.value numeric, a p-value
#'
#' @return character, an asterick associated with the given
#' p-value. < 0.001 yields `***`, < 0.01 yields `**``,
#'  < 0.05 yields `*`, and < 0.001 yields ``.
#' @export
#'
set_sig_level <- function(p.value) {
  dplyr::case_when(p.value < 0.001 ~ "***",
                   p.value < 0.01 ~ "**",
                   p.value < 0.05 ~ "*",
                   p.value < 0.001 ~ "")
}

#' Create a standard hue of colors for ggplot2
#'
#' @param n integer, number of colors to return for plotting.
#'
#' @return character, a hexcode set of colors of length n.
#' @export
#'
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length.out=n+1)
  grDevices::hcl(h=hues, l=65, c=100)[seq_len(n)]
}
