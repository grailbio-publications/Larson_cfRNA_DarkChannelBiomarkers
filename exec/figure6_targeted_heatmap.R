source('./config.R')

pilot_dcbs <- union(pilot_breast_dcbs, pilot_lung_dcbs) %>%
  union(heteroDE)

parse_arguments <- function() {

  parser <- argparse::ArgumentParser()

  parser$add_argument('--count_cfrna',
                      required = TRUE,
                      default = 'whole_transcriptome_cfrna_spliced_counts.tsv',
                      help = 'TSV file with whole transcriptome plasma gene counts')
  parser$add_argument('--anno_validation',
                      required = TRUE,
                      default = 'validation_hidden_plasma_ids.tsv',
                      help = '')
  parser$add_argument('--rpm_validation',
                      required = TRUE,
                      default = 'validation_plasma_rpms.tsv',
                      help = '')
  parser$add_argument('--count_validation',
                      required = TRUE,
                      default = 'validation_plasma_strict_counts.tsv',
                      help = '')

  args <- parser$parse_args()
  return(args)
}

main <- function() {

  args <- parse_arguments()

  # human protein atlas tissue specificity
  human_protein_atlas <- readr::read_tsv(system.file("extdata",
                                                     "human_protein_atlas_tissue_specificity.tsv",
                                                     package = "cellfreetranscriptome"))

  validation_metadata <- readr::read_tsv(args$anno_validation) %>%
    # arrange by indication and staging
    dplyr::arrange(indication, stage)
  validation_strict_counts <- readr::read_tsv(args$count_validation)
  validation_rpms <- readr::read_tsv(args$rpm_validation)

  # get norm, brca, luca samples
  norm_spls <- validation_metadata %>%
    dplyr::filter(indication == 'NORM') %>%
    dplyr::pull(`Sample ID`)
  luca_spls <- validation_metadata %>%
    dplyr::filter(indication %in% c('LUCA', 'LUAD')) %>%
    dplyr::pull(`Sample ID`)
  brca_spls <- validation_metadata %>%
    dplyr::filter(indication == 'BRCA') %>%
    dplyr::pull(`Sample ID`)

  # collapsed_count data frames for DE
  # edger expects raw counts, no need to normalize to rpms
  # edger normalization accounts for library size normalization
  norm_counts <- validation_strict_counts %>%
    dplyr::select(gene_id, norm_spls)
  brca_counts <- validation_strict_counts %>%
    dplyr::select(gene_id, brca_spls)
  luca_counts <- validation_strict_counts %>%
    dplyr::select(gene_id, luca_spls)

  norm_luca_counts <- norm_counts %>% dplyr::full_join(luca_counts)
  norm_brca_counts <- norm_counts %>% dplyr::full_join(brca_counts)

  # get darknes index
  darkness_index <- get_darkness_index(args$count_cfrna)

  luca_labels <- c(rep('Norm', ncol(norm_counts)-1),
                   rep('Luca', ncol(luca_counts)-1))
  ## luca ranked genes
  de_norm_luca <- edger_two_group_de_te(norm_luca_counts,
                                        luca_labels,
                                        c(biscotti_pos_control_genes)) %>%
    dplyr::filter(gene_id %in% union(biscotti_pos_control_genes,
                                     pilot_dcbs))

  brca_labels <- c(rep('Norm', ncol(norm_counts)-1),
                   rep('Brca', ncol(brca_counts)-1))
  ## brca ranked genes
  de_norm_brca <- edger_two_group_de_te(norm_brca_counts,
                                        brca_labels,
                                        c(biscotti_pos_control_genes)) %>%
    dplyr::filter(gene_id %in% union(biscotti_pos_control_genes,
                                     pilot_dcbs))
  # allows for <1 false positive out of 36 genes tested
  # 1 / n_genes_tested
  fdrThresh <- 0.01

  p <- ggplot2::ggplot(de_norm_luca,
                       ggplot2::aes(x = reorder(gene_id, -FDR), y = FDR, fill = logFC)) +
    ggplot2::geom_bar(stat = 'identity') + ggplot2::theme_bw() + ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = fdrThresh, linetype = 'dashed') +
    ggplot2::labs(x = '', title = 'Lung cancer DE ranked genes') +
    ggplot2::scale_y_log10() +
    ggplot2::theme(text = ggplot2::element_text(size = 12))
  p
  ggplot2::ggsave('figure6_te_de_luca_ranked_genes.pdf', device = 'pdf',
                  width = 5, height = 6, dpi = 150)

  de_luca_genes <- de_norm_luca %>%
    dplyr::filter(FDR < fdrThresh) %>%
    dplyr::arrange(FDR)

  # allows for <1 false positive out of 35 genes tested
  # 1 / n_genes_tested
  fdrThresh <- 0.01

  p <- ggplot2::ggplot(de_norm_brca,
                       ggplot2::aes(x = reorder(gene_id, -FDR), y = FDR, fill = logFC)) +
    ggplot2::geom_bar(stat = 'identity') + ggplot2::theme_bw() + ggplot2::coord_flip() +
    ggplot2::labs(x = '', title = 'Breast cancer DE ranked genes') +
    ggplot2::geom_hline(yintercept = fdrThresh, linetype = 'dashed') +
    ggplot2::scale_y_log10() +
    ggplot2::theme(text = ggplot2::element_text(size = 12))
  p
  ggplot2::ggsave('figure6_te_de_brca_ranked_genes.pdf', device = 'pdf',
                  width = 5, height = 6, dpi = 150)

  de_brca_genes <- de_norm_brca %>%
    dplyr::filter(FDR < fdrThresh) %>%
    dplyr::arrange(FDR)

  feature_anno <- data.frame(Gene = validation_rpms$gene_id) %>%
    dplyr::left_join(human_protein_atlas) %>%
    dplyr::select(Gene, `Gene description`, `RNA tissue category`, Tissue_specific) %>%
    dplyr::arrange(Tissue_specific) %>%
    dplyr::mutate(Tissue_specific = ifelse(is.na(Tissue_specific),
                                           'Non-specific',
                                           Tissue_specific),
                  Tissue_specific = as.factor(Tissue_specific)) %>%
    dplyr::select(Gene, Tissue_specific)

  dat <- validation_rpms %>%
    dplyr::filter(gene_id %in% pilot_dcbs) %>%
    dplyr::left_join(feature_anno, by = c('gene_id' = 'Gene')) %>%
    dplyr::arrange(Tissue_specific) %>%
    plyr::rename(c('gene_id' = 'Gene',
                   'Tissue_specific' = 'Tissue Specificity'))
  col_labels <- data.frame(sample = colnames(dat %>%
                                               dplyr::select(-c(Gene,
                                                                `Tissue Specificity`)))) %>%
    dplyr::mutate(`Cancer Type` = dplyr::case_when(sample %in% colnames(luca_counts) ~
                                                     'Lung',
                                                   sample %in% colnames(brca_counts) ~
                                                     'Breast',
                                                   sample %in% colnames(norm_counts) ~
                                                     'Non-cancer')) %>%
    dplyr::left_join(validation_metadata,
                     by = c('sample' = 'Sample ID')) %>%
    dplyr::mutate(stage = ifelse(stage == 'NA', NA, stage)) %>%
    dplyr::arrange(`Cancer Type`, stage) %>%
    dplyr::select(sample, stage, `Cancer Type`) %>%
    tibble::column_to_rownames('sample')
  my_colour <- list(stage = c(I = "grey93", II = "grey74",
                              III = "grey33", IV = "grey5"))

  row_labels <- dat %>%
    dplyr::select(`Tissue Specificity`, Gene) %>%
    dplyr::left_join(darkness_index, by = c('Gene' = 'gene')) %>%
    dplyr::arrange(`Tissue Specificity`, median, sd) %>%
    dplyr::select(sd, `Tissue Specificity`, Gene) %>%
    tibble::column_to_rownames('Gene')
  dat <- dat %>%
    dplyr::left_join(darkness_index, by = c('Gene' = 'gene')) %>%
    dplyr::arrange(`Tissue Specificity`, median, sd) %>%
    tibble::column_to_rownames('Gene') %>%
    dplyr::select(-c(`Tissue Specificity`, median, sd))
  dat <- dat[, rownames(col_labels)]

  pheatmap::pheatmap(t(scale(t(as.matrix(dat)), center = F)) %>% replace(is.na(.), 0),
                     scale = 'none',
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     color = colorRampPalette(RColorBrewer::brewer.pal(n = 7,
                                                                       name ="Purples"))(100),
                     annotation_row = row_labels,
                     annotation_col = col_labels,
                     annotation_colors = my_colour,
                     show_colnames = FALSE,
                     annotation_names_row = F,
                     annotation_name_row = T,
                     annotation_names_col = FALSE,
                     fontsize = 12,
                     gaps_row = c(4, 10),
                     gaps_col = c(35, 53),
                     annotation_legend = T,
                     border_color = NA,
                     width = 9, height = 6, filename = 'figure6_targeted_enriched_heatmap.pdf')

}

if (sys.nframe() == 0) {
  main()
}
