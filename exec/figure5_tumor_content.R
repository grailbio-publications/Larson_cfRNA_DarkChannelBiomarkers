source('./config.R')

parse_arguments <- function() {

  parser <- argparse::ArgumentParser()

  parser$add_argument('--count_cfrna',
                      required = TRUE,
                      default = 'whole_transcriptome_cfrna_spliced_counts.tsv',
                      help = 'TSV file with whole transcriptome plasma gene counts')
  parser$add_argument('--count_tissue',
                      required = TRUE,
                      default = 'whole_transcriptome_tissue_spliced_counts.tsv',
                      help = 'TSV file with whole transcriptome tissue gene counts')
  parser$add_argument('--tumor_fraction',
                      required = TRUE,
                      default = 'tumor_fraction_estimates.tsv',
                      help = 'TSV file with tumor fraction estimates for each participant')

  args <- parser$parse_args()
  return(args)
}

main <- function() {

  args <- parse_arguments()
  # get ccga plasma counts
  count_cfrna_df <- readr::read_tsv(args$count_cfrna) %>%
    dplyr::select(gene_id, dplyr::starts_with('BrCa'))
  sample_anno_cfrna <- data.frame(sample_id = colnames(count_cfrna_df)[-1]) %>%
    dplyr::mutate(class = dplyr::case_when(grepl('BrCa', sample_id) ~ 'Breast',
                                           grepl('LuCa', sample_id) ~ 'Lung',
                                           grepl('NC', sample_id) ~ 'Non-cancer'))

  # get ccga tissue counts
  count_tissue_df <- readr::read_tsv(args$count_tissue) %>%
    dplyr::select(gene_id, dplyr::starts_with('BrCa'))
  sample_anno_tissue <- data.frame(sample_id = colnames(count_tissue_df)[-1]) %>%
    dplyr::mutate(class = dplyr::case_when(grepl('BrCa', sample_id) ~ 'Breast',
                                           grepl('LuCa', sample_id) ~ 'Lung',
                                           grepl('NC', sample_id) ~ 'Non-cancer'))

  tumor_frac <- readr::read_tsv(args$tumor_fraction)

  # Create a sample list arranged by tumor fraction (highest to lowest)
  brca_sample_TF_order <- tumor_frac %>%
    dplyr::filter(primccat == "Breast") %>%
    tidyr::drop_na() %>%
    dplyr::arrange(tumor_fraction)

  # transpose data matrix
  count_cfrna_t <- genetable_to_matrix(count_cfrna_df, transpose=F)
  count_tissue_t <- genetable_to_matrix(count_tissue_df, transpose=F)
  # normalize read counts to library size
  norm_count_cfrna <- count_normalization(count_cfrna_t, sample_anno_cfrna)
  norm_count_tissue <- count_normalization(count_tissue_t, sample_anno_tissue)

  # SCGB2A2 plot
  sample_anno_breast <- sample_anno_cfrna %>%
    dplyr::filter(sample_id %in% brca_sample_TF_order$paper_id)

  gene <- "SCGB2A2"
  cfrna_expr <- as.numeric(norm_count_cfrna[gene, match(sample_anno_breast$sample_id,
                                                   colnames(norm_count_cfrna))])
  tissue_expr <- as.numeric(norm_count_tissue[gene, match(sample_anno_breast$sample_id,
                                                     colnames(norm_count_tissue))])
  plot.matrix <- cbind(sample_anno_breast, cfrna_expr, tissue_expr)
  colnames(plot.matrix)[ncol(plot.matrix) - 1] <- "rpm_cfrna"
  colnames(plot.matrix)[ncol(plot.matrix)] <- "rpm_tissue"
  plot.matrix <- plot.matrix %>%
    dplyr::left_join(tumor_frac,
                     by = c('sample_id' = 'paper_id')) %>%
    dplyr::mutate(tumor_content = tumor_fraction * rpm_tissue,
                  tumor_content_lcb = lcb * rpm_tissue,
                  tumor_content_ucb = ucb * rpm_tissue)%>%
    dplyr::mutate(
      detected = dplyr::case_when(
        .$rpm_cfrna > 0   ~ "Yes",
        .$rpm_cfrna <= 0  ~ "No",
        TRUE                      ~  "pass"
      )) %>%
    dplyr::arrange(tumor_content)

  color_group <- c("red", "gray60")
  names(color_group) <- c("Yes", "No")

  # sample order based on tumor fractions
  plot.matrix <- plot.matrix %>%
    dplyr::mutate(sample_id = factor(sample_id, levels = brca_sample_TF_order$paper_id))

  # SCGB2A2 expression plotted by patient IDs
  # as a function of tumor fraction (triangles) and tumor content (squares)
  p <-  ggplot2::ggplot(data=plot.matrix, ggplot2::aes(x=sample_id, y=tumor_content,
                                                       ymin=tumor_content_lcb,
                                                       ymax=tumor_content_ucb,
                                                       color = detected)) +
    ggplot2::geom_segment(ggplot2::aes(x=sample_id, xend=sample_id,
                                       y=tumor_fraction, yend=tumor_content), color="grey") +
    ggplot2::geom_point(ggplot2::aes(x=sample_id, y=tumor_fraction),
                        shape="triangle", size=2) +
    ggplot2::geom_point(ggplot2::aes(x=sample_id, y=tumor_content),
                        shape="square", size=2) +
    ggplot2::geom_hline(yintercept=1e-2) +
    ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
    ggplot2::xlab("patient id") +
    ggplot2::ylab("Tumor content") +
    ggplot2::scale_y_log10()+
    ggplot2::ggtitle(paste0("gene: ", gene))+
    ggplot2::scale_color_manual(values=color_group) +
    ggplot2::theme_bw(base_size = 6)+
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank())

  pdf(file = "./figure5_lollipop_SCGB2A2.pdf",
      width = 3.5, height =3.5,
      bg = "white")
  print(p)
  dev.off()

  # FABP7 plot
  gene <- "FABP7"
  cfrna_expr <- as.numeric(norm_count_cfrna[gene, match(sample_anno_breast$sample_id,
                                                   colnames(norm_count_cfrna))])
  tissue_expr <- as.numeric(norm_count_tissue[gene, match(sample_anno_breast$sample_id,
                                                     colnames(norm_count_tissue))])
  plot.matrix <- cbind(sample_anno_breast, cfrna_expr, tissue_expr)
  colnames(plot.matrix)[ncol(plot.matrix) - 1] <- "rpm_cfrna"
  colnames(plot.matrix)[ncol(plot.matrix)] <- "rpm_tissue"
  plot.matrix <- plot.matrix %>%
    dplyr::left_join(tumor_frac,
                     by = c('sample_id' = 'paper_id')) %>%
    dplyr::mutate(tumor_content = tumor_fraction * rpm_tissue,
                  tumor_content_lcb = lcb * rpm_tissue,
                  tumor_content_ucb = ucb * rpm_tissue)%>%
    dplyr::mutate(
      detected = dplyr::case_when(
        .$rpm_cfrna > 0   ~ "Yes",
        .$rpm_cfrna <= 0  ~ "No",
        TRUE                      ~  "pass"
      )) %>%
    dplyr::arrange(tumor_content)

  color_group <- c("blue", "gray60")
  names(color_group) <- c("Yes", "No")

  # sample order based on tumor fractions
  plot.matrix <- plot.matrix %>%
    dplyr::mutate(sample_id = factor(sample_id, levels = brca_sample_TF_order$paper_id))

  # as a function of tumor fraction (triangles) and tumor content (squares)
  p <-  ggplot2::ggplot(data=plot.matrix, ggplot2::aes(x=sample_id, y=tumor_content,
                                                       ymin=tumor_content_lcb,
                                                       ymax=tumor_content_ucb,
                                                       color = detected)) +
    ggplot2::geom_segment(ggplot2::aes(x=sample_id, xend=sample_id,
                                       y=tumor_fraction, yend=tumor_content), color="grey") +
    ggplot2::geom_point(ggplot2::aes(x=sample_id, y=tumor_fraction),
                        shape="triangle", size=2) +
    ggplot2::geom_point(ggplot2::aes(x=sample_id, y=tumor_content),
                        shape="square", size=2) +
    ggplot2::geom_hline(yintercept=1e-2) +
    ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
    ggplot2::xlab("patient id") +
    ggplot2::ylab("Tumor content") +
    ggplot2::scale_y_log10()+
    ggplot2::ggtitle(paste0("gene: ", gene))+
    ggplot2::scale_color_manual(values=color_group) +
    ggplot2::theme_bw(base_size = 6)+
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank())
  pdf(file = "./figure5_lollipop_FABP7.pdf",
      width = 3.5, height =3.5,
      bg = "white")
  print(p)
  dev.off()

}

if (sys.nframe() == 0) {
  main()
}
