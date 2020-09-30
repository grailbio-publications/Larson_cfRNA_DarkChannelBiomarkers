source('./config.R')

parse_arguments <- function() {

  parser <- argparse::ArgumentParser()

  parser$add_argument('--gtex_median',
                      required = TRUE,
                      default = 'gtex_rpm_median.RData',
                      help = 'RData file with GTEx median RPMs')
  parser$add_argument('--count_cfrna',
                      required = TRUE,
                      default = 'whole_transcriptome_cfrna_spliced_counts.tsv',
                      help = 'TSV file with whole transcriptome plasma gene counts')
  parser$add_argument('--count_tissue',
                      required = TRUE,
                      default = 'whole_transcriptome_tissue_spliced_counts.tsv',
                      help = 'TSV file with whole transcriptome tissue gene counts')

  args <- parser$parse_args()
  return(args)
}

main <- function() {

  args <- parse_arguments()

  # load preprocessed GTEx samples
  load(args$gtex_median)
  gtex_symbol_rpm_median <- as.data.frame(gtex_rpm_median) %>%
    tibble::rownames_to_column(., "gene_id") %>%
    convert_gene_id()

  # compute tissue specificity score
  gtex_tss <- compute_TSS(gtex_rpm_median) %>%
    convert_gene_id()

  ## Load and preprocess the CCGA cfRNA data
  count_cfrna_df <- readr::read_tsv(args$count_cfrna)
  count_tissue_df <- readr::read_tsv(args$count_tissue)

  # Reads count normalization, cfrna
  rpm_cfrna_symbol <- get_rpm(count_cfrna_df) %>%
    # log2 RPM
    dplyr::mutate_if(is.numeric, function(x) log2(x+1))
  # Reads count normalization, tissue
  rpm_tissue_symbol <- get_rpm(count_tissue_df) %>%
    # log2 RPM
    dplyr::mutate_if(is.numeric, function(x) log2(x+1))

  # Only consider the mutual genes between CCGA and GTEx datasets (n=17071)
  mutual_genes <- intersect(gtex_symbol_rpm_median$gene_id,
                            rpm_cfrna_symbol$gene_id)

  TSS_matrix_GTEx_mutual <- gtex_tss %>%
    dplyr::filter(gene_id %in% mutual_genes)

  gtex_rpm_median_log <- gtex_symbol_rpm_median %>%
    dplyr::filter(gene_id %in% mutual_genes) %>%
    dplyr::mutate_if(is.numeric, function(x) log2(x+1)) %>%
    tibble::column_to_rownames("gene_id")

  # Prepare the tissue signatures
  ## Set hyperparameters
  TSS_lowerBound <- 0.40
  TSS_higherBound <- 1
  # the number signature genes per tissue
  gene_per_tissue <- 20

  tissue_base <- unique(TSS_matrix_GTEx_mutual$Top_tissue)

  ## select the signature genes for each tissue.
  TSS_matrix_signature <- c()
  for (tissue in tissue_base) {
    TSS_matrix_signature_singleTissue <- TSS_matrix_GTEx_mutual %>%
      dplyr::filter(Top_tissue == tissue &
                      Tissue_expr > 50 &
                      Tissue_expr < 10000 &
                      Expr_weight > TSS_lowerBound) %>%
      dplyr::arrange(desc(Expr_weight)) %>%
      dplyr::slice(seq_len(gene_per_tissue))
    TSS_matrix_signature <- rbind(TSS_matrix_signature,
                                  TSS_matrix_signature_singleTissue)
  }

  genes_vec <- TSS_matrix_signature$gene_id
  signatures <- as.data.frame(gtex_rpm_median_log[match(genes_vec,
                                                        rownames(gtex_rpm_median_log)),
                                                  match(tissue_base,
                                                        colnames(gtex_rpm_median_log))])
  n_genes <- nrow(signatures)

  # Tissue deconvolution of bulk RNAseq collected from non-cancer plasma
  # Quadratic programming on the CCGA cfrna samples
  sample_selected <- colnames(rpm_cfrna_symbol[-1])
  dataset_DeconRNASeq <- as.data.frame(rpm_cfrna_symbol[match(genes_vec,
                                                              rpm_cfrna_symbol$gene_id),
                                                        sample_selected])

  # normalize base and tissue matrices
  signatures_norm <- normalize_matrix(as.matrix(signatures))
  dataset_DeconRNASeq_norm <- normalize_matrix(as.matrix(dataset_DeconRNASeq))

  # deconvolute and merge back to annotation
  fraction_predicted_cfrna <- tissueDeconResidual(signatures_norm,
                                                  dataset_DeconRNASeq_norm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column('sample_id') %>%
    dplyr::mutate(type = dplyr::case_when(grepl('BrCa', sample_id) ~ 'Breast',
                                           grepl('LuCa', sample_id) ~ 'Lung',
                                           grepl('NC', sample_id) ~ 'Non-cancer'))

  cols <- colnames(fraction_predicted_cfrna)[2:31]
  goodness_of_fit <- apply(fraction_predicted_cfrna, 1, function(x) {
    spl <- unlist(x['sample_id'])
    vec <- as.numeric(matrix(unlist(x[cols]), ncol=1, nrow=30))
    corr <- cor(matrix(unlist(signatures), ncol = 30, nrow=n_genes) %*% vec,
                matrix(unlist(dataset_DeconRNASeq[spl])))
    data.frame(sample_id = spl,
               type = unlist(x['type']),
               pearsons = corr)
  }) %>%
    do.call('rbind', .)

  # order samples by decreasing blood fraction
  spl_order <- unlist(fraction_predicted_cfrna$sample_id[order(-fraction_predicted_cfrna$Blood)])
  dat <- fraction_predicted_cfrna %>%
    tidyr::gather(tissue_type, percentage, -type, -sample_id) %>%
    dplyr::filter(type %in% c('Non-cancer', 'Breast', 'Lung')) %>%
    tidyr::spread(tissue_type, percentage) %>%
    tidyr::gather(tissue_type, percentage, -type, -sample_id)

  nc_dat <- dat %>%
    dplyr::filter(type == 'Non-cancer')  %>%
    dplyr::arrange(tissue_type, -percentage)

  # median blood content across samples
  blood_median <- quantile(nc_dat[which(nc_dat$tissue_type == 'Blood'), ]$percentage, 0.5)

  ##Stacked bar chart for non-cancer samples
  p_nc_plasma <- nc_dat %>%
    dplyr::mutate(tissue_type=ifelse(tissue_type %in% c("Blood", "Spleen",
                                                         "Liver"),
                                     tissue_type,
                                     "Remainder"),
                  sample_id = factor(sample_id, levels = spl_order)) %>%
    ggplot2::ggplot(ggplot2::aes(x=sample_id, y=percentage,
                                 fill=reorder(tissue_type, -percentage),
                                 label=tissue_type))+
    ggplot2::geom_bar(stat="identity",
                      position = ggplot2::position_stack(reverse = TRUE),
                      alpha = 0.7) +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 20),
                   legend.position = c(0.15, 0.2)) +
    ggplot2::geom_hline(yintercept = blood_median, linetype = 'dashed') +
    ggplot2::labs(y= "Tissue fraction", x = "") +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis(name = '')) +
    ggplot2::scale_fill_manual(values = c("darkorchid1", "orange1",
                                          "darkslategray1", "chartreuse4"))

  # Set color scheme
  type <- c("BrCa", "LuCa", "Non-cancer")
  num_of_type <- length(type) - 1
  color_scheme <- c(gg_color_hue(num_of_type), "gray60")
  names(color_scheme) <- type


  p_lung_fraction_cancer_plasma <- fraction_predicted_cfrna %>%
    tidyr::gather(key, value, -type, -sample_id) %>%
    dplyr::filter(type %in% c('Non-cancer', 'Lung', 'Breast'),
                  key == 'Lung') %>%
    dplyr::mutate(
      type = dplyr::case_when(type == 'Non-cancer' ~ 'Non-cancer',
                               type == 'Lung' ~ 'LuCa',
                               type == 'Breast' ~ 'BrCa'),
      type = factor(type, levels = c('Non-cancer', 'BrCa', 'LuCa'))) %>%
    ggplot2::ggplot(ggplot2::aes(x = type, y = value, colour = type)) +
    ggplot2::geom_boxplot(outlier.shape = NA, notch = F) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.25,
                                                                 jitter.height= 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 20),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::labs(x= '', y = 'Lung fraction') +
    ggplot2::scale_color_manual(values=color_scheme)

  p_lung_fraction_cancer_plasma

  p_breast_fraction_cancer_plasma <- fraction_predicted_cfrna %>%
    tidyr::gather(key, value, -type, -sample_id) %>%
    dplyr::filter(type %in% c('Non-cancer', 'Lung', 'Breast'),
                  key == 'Breast') %>%
    dplyr::mutate(
      type = dplyr::case_when(type == 'Non-cancer' ~ 'Non-cancer',
                               type == 'Lung' ~ 'LuCa',
                               type == 'Breast' ~ 'BrCa'),
      type = factor(type, levels = c('Non-cancer', 'BrCa', 'LuCa'))) %>%
    ggplot2::ggplot(ggplot2::aes(x = type, y = value, colour = type)) +
    ggplot2::geom_boxplot(outlier.shape = NA, notch = F) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.25,
                                                                 jitter.height= 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 20),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::labs(x= '', y = 'Breast fraction') +
    ggplot2::scale_color_manual(values=color_scheme)

  # Tissue deconvolution of bulk RNAseq collected from tumor tissue
  # subset tissue expression matrix by genes of interest
  dataset_DeconRNASeq_tissue <- rpm_tissue_symbol[match(genes_vec,
                                                        rpm_tissue_symbol$gene_id), ]
  dataset_DeconRNASeq_tissue <- dataset_DeconRNASeq_tissue[-1]

  # normalize base and tissue matrices
  signatures_norm <- normalize_matrix(as.matrix(signatures))
  dataset_DeconRNASeq_norm <- normalize_matrix(as.matrix(dataset_DeconRNASeq_tissue))

  # deconvolute and merge back to annotation
  fraction_predicted_tissue <- tissueDeconResidual(signatures_norm,
                                                   dataset_DeconRNASeq_norm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column('sample_id') %>%
    dplyr::mutate(type = dplyr::case_when(grepl('BrCa', sample_id) ~ 'Breast',
                                           grepl('LuCa', sample_id) ~ 'Lung'))

  p_breast_fraction_cancer_tissue <- fraction_predicted_tissue %>%
    tidyr::gather(key, value, -type, -sample_id) %>%
    dplyr::filter(type %in% c('Lung', 'Breast'),
                  key == 'Breast') %>%
    dplyr::mutate(
      type = dplyr::case_when(type == 'Lung' ~ 'LuCa',
                               type == 'Breast' ~ 'BrCa'),
      type = factor(type, levels = c('BrCa', 'LuCa'))) %>%
    ggplot2::ggplot(ggplot2::aes(x = type, y = value, colour = type)) +
    ggplot2::geom_boxplot(outlier.shape = NA, notch = F) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.25,
                                                                 jitter.height= 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 20),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::labs(x= '', y = 'Breast fraction') +
    ggplot2::scale_color_manual(values=color_scheme)

  p_lung_fraction_cancer_tissue <- fraction_predicted_tissue %>%
    tidyr::gather(key, value, -type, -sample_id) %>%
    dplyr::filter(type %in% c('Lung', 'Breast'),
                  key == 'Lung') %>%
    dplyr::mutate(
      type = dplyr::case_when(type == 'Lung' ~ 'LuCa',
                               type == 'Breast' ~ 'BrCa'),
      type = factor(type, levels = c('BrCa', 'LuCa'))) %>%
    ggplot2::ggplot(ggplot2::aes(x = type, y = value, colour = type)) +
    ggplot2::geom_boxplot(outlier.shape = NA, notch = F) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.25,
                                                                 jitter.height= 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 20),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::labs(x= '', y = 'Lung fraction') +
    ggplot2::scale_color_manual(values=color_scheme)

  # Arrange
  title_noncancer <- cowplot::ggdraw() + cowplot::draw_label("Non-cancer plasma deconvolution",
                                                             size = 30, y = 0.4)
  title_tissue <- cowplot::ggdraw() + cowplot::draw_label("Tumor tissue", size = 30, y = 0.4)
  title_plasma <- cowplot::ggdraw() + cowplot::draw_label("Plasma", size = 30, y = 0.4)

  col1 <- cowplot::plot_grid(title_noncancer, p_nc_plasma,
                             ncol=1, rel_heights=c(0.1, 1))
  col2 <- cowplot::plot_grid(title_tissue,
                             p_lung_fraction_cancer_tissue +
                               ggplot2::theme(legend.position = 'none'),
                             p_breast_fraction_cancer_tissue +
                               ggplot2::theme(legend.position = 'none'),
                             ncol=1, rel_heights=c(0.1, 0.5, 0.5))
  col3 <- cowplot::plot_grid(title_plasma,
                             p_lung_fraction_cancer_plasma +
                               ggplot2::theme(legend.position = 'none'),
                             p_breast_fraction_cancer_plasma +
                               ggplot2::theme(legend.position = 'none'),
                             ncol=1, rel_heights=c(0.1, 0.5, 0.5))

  cowplot::plot_grid(col1,
                     col2,
                     col3,
                     ncol = 3,
                     rel_widths =c(1, 0.4, 0.6),
                     labels = c('A', 'B', 'C'),
                     label_size = 30,
                     align = 'hv')

  ggplot2::ggsave('figure2_deconvolution_final.pdf', device = 'pdf',
                  width = 25, height = 10, dpi = 300)


}

if (sys.nframe() == 0) {
  main()
}
