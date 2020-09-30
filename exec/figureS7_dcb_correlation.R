source('./config.R')

pilot_dcbs <- union(pilot_breast_dcbs, pilot_lung_dcbs) %>%
  union(heteroDE)

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
  parser$add_argument('--anno_cfrna',
                      required = TRUE,
                      default = 'whole_transcriptome_hidden_plasma_ids.tsv',
                      help = 'TSV file with whole transcriptome plasma annotations')
  parser$add_argument('--anno_tissue',
                      required = TRUE,
                      default = 'whole_transcriptome_hidden_tissue_ids.tsv',
                      help = 'TSV file with whole transcriptome tissue annotations')
  parser$add_argument('--tcga_luad',
                      required = TRUE,
                      default = 'luad_tumor_counts.tsv',
                      help = '')
  parser$add_argument('--tcga_hrpos_her2pos',
                      required = TRUE,
                      default = 'brca_hrpos_her2pos_tumor_counts.tsv',
                      help = '')
  parser$add_argument('--tcga_hrpos_her2neg',
                      required = TRUE,
                      default = 'brca_hrpos_her2neg_tumor_counts.tsv',
                      help = '')
  parser$add_argument('--tcga_hrneg_her2pos',
                      required = TRUE,
                      default = 'brca_hrneg_her2pos_tumor_counts.tsv',
                      help = '')
  parser$add_argument('--tcga_tnbc',
                      required = TRUE,
                      default = 'brca_tnbc_tumor_counts.tsv',
                      help = '')


  args <- parser$parse_args()
  return(args)
}

main <- function() {

  args <- parse_arguments()
  # get darknes index
  darkness_index <- get_darkness_index(args$count_cfrna)
  # get dark channels
  dark_channels <- darkness_index %>%
    dplyr::filter(median < 1, sd < 0.1) %>%
    dplyr::select(gene) %>%
    unlist()

  # get ccga tissue counts for lung and breast
  count_tissue_df <- readr::read_tsv(args$count_tissue)
  ccga_breast_tissue_rpm <- count_tissue_df %>%
    dplyr::select(gene_id, dplyr::starts_with('BrCa')) %>%
    get_rpm()
  ccga_lung_tissue_rpm <- count_tissue_df %>%
    dplyr::select(gene_id, dplyr::starts_with('LuCa')) %>%
    get_rpm()

  # get ccga plasma counts for lung and breast
  count_plasma_df <- readr::read_tsv(args$count_cfrna)
  ccga_breast_plasma_rpm <- count_plasma_df %>%
    dplyr::select(gene_id, dplyr::starts_with('BrCa')) %>%
    get_rpm()
  ccga_lung_plasma_rpm <- count_plasma_df %>%
    dplyr::select(gene_id, dplyr::starts_with('LuCa')) %>%
    get_rpm()

  # get tcga gene counts
  tcga_lung_tissue_rpm <- readr::read_tsv(args$tcga_luad) %>%
    get_rpm()
  brca_hrneg_her2pos_rpm <- readr::read_tsv(args$tcga_hrneg_her2pos) %>%
    get_rpm()
  brca_hrpos_her2neg_rpm <- readr::read_tsv(args$tcga_hrpos_her2neg) %>%
    get_rpm()
  brca_hrpos_her2pos_rpm <- readr::read_tsv(args$tcga_hrpos_her2pos) %>%
    get_rpm()
  brca_tnbc_rpm <- readr::read_tsv(args$tcga_tnbc) %>%
    get_rpm()
  tcga_breast_tissue_rpm <- cbind(brca_hrneg_her2pos_rpm,
                                  brca_hrpos_her2neg_rpm) %>%
    cbind(brca_hrpos_her2pos_rpm) %>%
    cbind(brca_tnbc_rpm)

  # summarise gene counts
  ccga_lung_tissue_summary <- summarise_rpms(ccga_lung_tissue_rpm$gene_id,
                                             ccga_lung_tissue_rpm %>%
                                               dplyr::select(-gene_id) %>%
                                               as.matrix())
  ccga_breast_tissue_summary <- summarise_rpms(ccga_breast_tissue_rpm$gene_id,
                                               ccga_breast_tissue_rpm %>%
                                                 dplyr::select(-gene_id) %>%
                                                 as.matrix())
  ccga_lung_plasma_summary <- summarise_rpms(ccga_lung_plasma_rpm$gene_id,
                                             ccga_lung_plasma_rpm %>%
                                               dplyr::select(-gene_id) %>%
                                               as.matrix())
  ccga_breast_plasma_summary <- summarise_rpms(ccga_breast_plasma_rpm$gene_id,
                                               ccga_breast_plasma_rpm %>%
                                                 dplyr::select(-gene_id) %>%
                                                 as.matrix())
  tcga_lung_summary <- summarise_rpms(tcga_lung_tissue_rpm$gene_id,
                                      tcga_lung_tissue_rpm %>%
                                        dplyr::select(-gene_id) %>%
                                        as.matrix())
  tcga_breast_summary <- summarise_rpms(tcga_breast_tissue_rpm$gene_id,
                                        tcga_breast_tissue_rpm %>%
                                          dplyr::select(-gene_id) %>%
                                          as.matrix())

  # Dark channel expression in CCGA plasma correlated with CCGA tumor tissue expression
  # Genes which have mean plasma or tissue expression of 0 are assigned pseudocounts of 1e-4
  spear_p1 <- dplyr::full_join(ccga_breast_plasma_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               ccga_breast_tissue_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               by = 'gene_id', suffix = c('_plasma', '_tissue')) %>%
    dplyr::summarise(cor = cor(mean_rpm_tissue, mean_rpm_plasma, method = 'spearman',
                               use = 'complete.obs')) %>%
    unlist()

  p1 <- dplyr::full_join(ccga_breast_plasma_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         ccga_breast_tissue_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         by = 'gene_id', suffix = c('_plasma', '_tissue')) %>%
    dplyr::mutate(mean_rpm_tissue = ifelse(mean_rpm_tissue ==0,
                                           1e-4,
                                           mean_rpm_tissue),
                  mean_rpm_plasma = ifelse(mean_rpm_plasma ==0,
                                           1e-4,
                                           mean_rpm_plasma),
                  gene_id = as.character(gene_id)) %>%
    ggplot2::ggplot(ggplot2::aes(x = mean_rpm_tissue, y = mean_rpm_plasma,
                                 color=ifelse(gene_id %in% pilot_breast_dcbs,
                                              'red', 'black'),
                                 label = gene_id)) +
    ggplot2::geom_point(alpha=0.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size =12),
                   legend.position = 'none') +
    ggplot2::scale_color_identity(guide = 'legend',
                                  name = 'pilot breast DCB',
                                  labels = c('FALSE', 'TRUE')) +
    ggrepel::geom_label_repel(ggplot2::aes(label=ifelse((gene_id %in% pilot_breast_dcbs),
                                                        gene_id,
                                                        '')),
                              show.legend = F) +
    ggplot2::labs(title = "Breast cancer",
                  y = 'CCGA plasma mean RPM',
                  x = 'CCGA tissue mean RPM') +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()


  spear_p2 <- dplyr::full_join(ccga_lung_plasma_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               ccga_lung_tissue_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               by = 'gene_id', suffix = c('_plasma', '_tissue')) %>%
    dplyr::summarise(cor = cor(mean_rpm_tissue, mean_rpm_plasma,
                               method = 'spearman', use = 'complete.obs')) %>%
    unlist()

  p2 <- dplyr::full_join(ccga_lung_plasma_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         ccga_lung_tissue_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         by = 'gene_id', suffix = c('_plasma', '_tissue')) %>%
    dplyr::mutate(mean_rpm_tissue = ifelse(mean_rpm_tissue ==0,
                                           5e-4,
                                           mean_rpm_tissue),
                  mean_rpm_plasma = ifelse(mean_rpm_plasma ==0,
                                           5e-4,
                                           mean_rpm_plasma),
                  gene_id = as.character(gene_id)) %>%
    ggplot2::ggplot(ggplot2::aes(x = mean_rpm_tissue, y = mean_rpm_plasma,
                                 color=ifelse(gene_id %in% pilot_lung_dcbs,
                                              'red', 'black'),
                                 label = as.character(gene_id))) +
    ggplot2::geom_point(alpha=0.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size =12),
                   legend.position = 'none') +
    ggplot2::scale_color_identity(guide = 'legend',
                                  name = 'pilot lung DCB',
                                  labels = c('FALSE', 'TRUE')) +
    ggrepel::geom_label_repel(ggplot2::aes(label=ifelse((gene_id %in% pilot_lung_dcbs),
                                                        as.character(gene_id),
                                                        '')),
                              show.legend = F) +
    ggplot2::labs(title = "Lung cancer",
                  y = 'CCGA plasma mean RPM',
                  x = 'CCGA tissue mean RPM') +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  # Dark channel expression in CCGA tumor tissue correlated with TCGA tumor tissue expression
  spear_p3 <- dplyr::full_join(ccga_breast_tissue_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               tcga_breast_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               by = 'gene_id', suffix = c('_ccga', '_tcga')) %>%
    dplyr::summarise(cor = cor(mean_rpm_tcga, mean_rpm_ccga,
                               method = 'spearman', use = 'complete.obs')) %>%
    unlist()

  p3 <- dplyr::full_join(ccga_breast_tissue_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         tcga_breast_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         by = 'gene_id',
                         suffix = c('_ccga', '_tcga')) %>%
    dplyr::mutate(mean_rpm_ccga = ifelse(mean_rpm_ccga ==0,
                                         1e-4,
                                         mean_rpm_ccga),
                  mean_rpm_tcga = ifelse(mean_rpm_tcga ==0,
                                         1e-5,
                                         mean_rpm_tcga),
                  gene_id = as.character(gene_id)) %>%
    ggplot2::ggplot(ggplot2::aes(x = mean_rpm_tcga, y = mean_rpm_ccga,
                                 color=ifelse(gene_id %in% pilot_breast_dcbs,
                                              'red', 'black'),
                                 label = as.character(gene_id))) +
    ggplot2::geom_point(alpha=0.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size =12),
                   legend.position = 'none') +
    ggplot2::scale_color_identity(guide = 'legend',
                                  name = 'pilot breast DCB',
                                  labels = c('FALSE', 'TRUE')) +
    ggrepel::geom_label_repel(ggplot2::aes(label=ifelse((gene_id %in% pilot_breast_dcbs),
                                                        as.character(gene_id), '')),
                              show.legend = F) +
    ggplot2::labs(title = paste0("Spearman's correlation coefficient: ", round(spear_p3, 3)),
                  y = 'CCGA tissue mean RPM',
                  x = 'TCGA tissue mean RPM') +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  spear_p4 <- dplyr::full_join(ccga_lung_tissue_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               tcga_lung_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               by = 'gene_id', suffix = c('_ccga', '_tcga')) %>%
    dplyr::summarise(cor = cor(mean_rpm_tcga, mean_rpm_ccga,
                               method = 'spearman', use = 'complete.obs')) %>%
    unlist()

  p4 <- dplyr::full_join(ccga_lung_tissue_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         tcga_lung_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         by = 'gene_id', suffix = c('_ccga', '_tcga')) %>%
    dplyr::mutate(mean_rpm_ccga = ifelse(mean_rpm_ccga ==0, 1e-4, mean_rpm_ccga),
                  mean_rpm_tcga = ifelse(mean_rpm_tcga ==0, 1e-5, mean_rpm_tcga),
                  gene_id = as.character(gene_id)) %>%
    ggplot2::ggplot(ggplot2::aes(x = mean_rpm_tcga, y = mean_rpm_ccga,
                                 color=ifelse(gene_id %in% pilot_lung_dcbs,
                                              'red', 'black'),
                                 label = as.character(gene_id))) +
    ggplot2::geom_point(alpha=0.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size =12),
                   legend.position = 'none') +
    ggplot2::scale_color_identity(guide = 'legend',
                                  name = 'pilot lung DCB',
                                  labels = c('FALSE', 'TRUE')) +
    ggrepel::geom_label_repel(ggplot2::aes(label=ifelse((gene_id %in% pilot_lung_dcbs),
                                                        as.character(gene_id), '')),
                              show.legend = F) +
    ggplot2::labs(title = paste0("Spearman's correlation coefficient: ", round(spear_p4, 3)),
                  y = 'CCGA tissue mean RPM',
                  x = 'TCGA tissue mean RPM') +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  # Dark channel expression in CCGA plasma correlated with TCGA tumor tissue expression
  spear_p5 <- dplyr::full_join(ccga_breast_plasma_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               tcga_breast_summary %>%
                                 dplyr::filter(gene_id %in% dark_channels),
                               by = 'gene_id', suffix = c('_ccga_plasma', '_tcga_tissue')) %>%
    dplyr::summarise(cor = cor(mean_rpm_tcga_tissue, mean_rpm_ccga_plasma,
                               method = 'spearman', use = 'complete.obs')) %>%
    unlist()

  p5 <- dplyr::full_join(ccga_breast_plasma_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         tcga_breast_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         by = 'gene_id',
                         suffix = c('_ccga_plasma', '_tcga_tissue')) %>%
    dplyr::mutate(mean_rpm_ccga_plasma = ifelse(mean_rpm_ccga_plasma ==0,
                                                1e-4,
                                                mean_rpm_ccga_plasma),
                  mean_rpm_tcga_tissue = ifelse(mean_rpm_tcga_tissue ==0,
                                                1e-5,
                                                mean_rpm_tcga_tissue)) %>%
    ggplot2::ggplot(ggplot2::aes(x = mean_rpm_tcga_tissue, y = mean_rpm_ccga_plasma,
                                 color=ifelse(gene_id %in% pilot_breast_dcbs,
                                              'red', 'black'),
                                 label = as.character(gene_id))) +
    ggplot2::geom_point(alpha=0.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size =12),
                   legend.position = 'none') +
    ggplot2::scale_color_identity(guide = 'legend',
                                  name = 'pilot breast DCB',
                                  labels = c('FALSE', 'TRUE')) +
    ggrepel::geom_label_repel(ggplot2::aes(label=ifelse((gene_id %in% pilot_breast_dcbs),
                                                        as.character(gene_id),
                                                        '')),
                              show.legend = F) +
    ggplot2::labs(title = paste0("Spearman's correlation coefficient: ", round(spear_p5, 3)),
                  y = 'CCGA plasma mean RPM',
                  x = 'TCGA tissue mean RPM') +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  spear_p6 <- dplyr::full_join(
    ccga_lung_plasma_summary %>%
      dplyr::filter(gene_id %in% dark_channels),
    tcga_lung_summary %>%
      dplyr::filter(gene_id %in% dark_channels),
    by = 'gene_id',
    suffix = c('_ccga_plasma', '_tcga_tissue')
  ) %>%
    dplyr::summarise(
      cor = cor(
        mean_rpm_tcga_tissue,
        mean_rpm_ccga_plasma,
        method = 'spearman',
        use = 'complete.obs'
      )
    ) %>%
    unlist()

  p6 <- dplyr::full_join(ccga_lung_plasma_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         tcga_lung_summary %>%
                           dplyr::filter(gene_id %in% dark_channels),
                         by = 'gene_id', suffix = c('_ccga_plasma', '_tcga_tissue')) %>%
    dplyr::mutate(mean_rpm_ccga_plasma = ifelse(mean_rpm_ccga_plasma ==0,
                                                1e-3,
                                                mean_rpm_ccga_plasma),
                  mean_rpm_tcga_tissue = ifelse(mean_rpm_tcga_tissue ==0,
                                                1e-5,
                                                mean_rpm_tcga_tissue)) %>%
    ggplot2::ggplot(ggplot2::aes(x = mean_rpm_tcga_tissue, y = mean_rpm_ccga_plasma,
                                 color=ifelse(gene_id %in% pilot_lung_dcbs,
                                              'red', 'black'),
                                 label = as.character(gene_id))) +
    ggplot2::geom_point(alpha=0.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size =12),
                   legend.position = 'none') +
    ggplot2::scale_color_identity(guide = 'legend',
                                  name = 'pilot lung DCB',
                                  labels = c('FALSE', 'TRUE')) +
    ggrepel::geom_label_repel(ggplot2::aes(label=ifelse((gene_id %in% pilot_lung_dcbs),
                                                        as.character(gene_id), '')),
                              show.legend = F) +
    ggplot2::labs(title = paste0("Spearman's correlation coefficient: ", round(spear_p6, 3)),
                  y = 'CCGA plasma mean RPM',
                  x = 'TCGA tissue mean RPM') +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  p_figs <- cowplot::plot_grid(p1, p2,
                               ncol=2, align = 'hv',
                               labels = c('AUTO'))

  p_figs
  ggplot2::ggsave('figureS7_tissue_plasma_correlation.pdf', device = 'pdf',
                  width = 14, height = 7, dpi = 300)

}

if (sys.nframe() == 0) {
  main()
}
