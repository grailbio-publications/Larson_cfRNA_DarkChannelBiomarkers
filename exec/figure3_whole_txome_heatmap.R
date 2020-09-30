source('./config.R')

pilot_dcbs <- union(pilot_breast_dcbs, pilot_lung_dcbs)

parse_arguments <- function() {

  parser <- argparse::ArgumentParser()

  parser$add_argument('--count_cfrna',
                      required = TRUE,
                      help = 'TSV file with whole transcriptome plasma gene counts')

  args <- parser$parse_args()
  return(args)
}

main <- function() {

  args <- parse_arguments()

  # WT heatmap
  # Plot the heatmap of RPMs from strict counts of whole tanscriptome pilot
  # Genes are ordered by darkness

  # human protein atlas tissue specificity
  human_protein_atlas <- readr::read_tsv(system.file("extdata",
                                                     "human_protein_atlas_tissue_specificity.tsv",
                                     package = "cellfreetranscriptome"))

  # get ccga plasma counts
  count_cfrna_df <- readr::read_tsv(args$count_cfrna)
  sample_anno_cfrna <- data.frame(sample_id = colnames(count_cfrna_df)[-1]) %>%
    dplyr::mutate(class = dplyr::case_when(grepl('BrCa', sample_id) ~ 'Breast',
                                           grepl('LuCa', sample_id) ~ 'Lung',
                                           grepl('NC', sample_id) ~ 'Non-cancer'))

  # sample IDs
  nc_spls <- sample_anno_cfrna %>%
    dplyr::filter(class == 'Non-cancer') %>%
    dplyr::pull(sample_id)
  breast_spls <- sample_anno_cfrna %>%
    dplyr::filter(class == 'Breast') %>%
    dplyr::pull(sample_id)
  lung_spls <- sample_anno_cfrna %>%
    dplyr::filter(class == 'Lung') %>%
    dplyr::pull(sample_id)

  # sample counts
  noncancer_counts <- count_cfrna_df %>%
    dplyr::select(gene_id, dplyr::starts_with('NC'))
  breast_counts <- count_cfrna_df %>%
    dplyr::select(gene_id, dplyr::starts_with('BrCa'))
  lung_counts <- count_cfrna_df %>%
    dplyr::select(gene_id, dplyr::starts_with('LuCa'))

  # sample RPMs
  rpm_nc <- get_rpm(noncancer_counts)
  rpm_breast <- get_rpm(breast_counts)
  rpm_lung <- get_rpm(lung_counts)

  # get darknes index
  darkness_index <- get_darkness_index(args$count_cfrna)

  rpm_df_wt <- rpm_breast %>%
    dplyr::full_join(rpm_lung, by = 'gene_id') %>%
    dplyr::full_join(rpm_nc, by = 'gene_id') %>%
    replace(is.na(.), 0)


  feature_anno <- data.frame(Gene = rpm_df_wt$gene_id) %>%
    dplyr::left_join(human_protein_atlas) %>%
    dplyr::select(Gene, `Gene description`, `RNA tissue category`, Tissue_specific) %>%
    dplyr::arrange(Tissue_specific) %>%
    dplyr::mutate(Tissue_specific = ifelse(is.na(Tissue_specific),
                                           'Non-specific',
                                           Tissue_specific),
                  Tissue_specific = as.factor(Tissue_specific)) %>%
    dplyr::select(Gene, Tissue_specific)

  dat <- rpm_df_wt %>%
    dplyr::filter(gene_id %in% pilot_dcbs) %>%
    dplyr::left_join(feature_anno, by = c('gene_id' = 'Gene')) %>%
    dplyr::arrange(Tissue_specific) %>%
    plyr::rename(c('gene_id' = 'Gene',
                   'Tissue_specific' = 'Tissue Specificity'))

  col_labels <- data.frame(sample = colnames(dat %>% dplyr::select(-c(Gene,
                                                                      `Tissue Specificity`)))) %>%
    dplyr::mutate(`Cancer Type` = dplyr::case_when(sample %in% colnames(rpm_lung) ~ 'Lung',
                                                   sample %in% colnames(rpm_breast) ~ 'Breast',
                                                   sample %in% colnames(rpm_nc) ~ 'Non-cancer')) %>%
    tibble::column_to_rownames('sample')
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

  pheatmap::pheatmap(t(scale(t(as.matrix(dat)), center = F)),
                     scale = 'none',
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     color = colorRampPalette(RColorBrewer::brewer.pal(n = 7,
                                                                       name ="Purples"))(100),
                     annotation_row = row_labels,
                     annotation_col = col_labels,
                     show_colnames = FALSE,
                     annotation_names_row = F,
                     annotation_name_row = T,
                     annotation_names_col = FALSE,
                     fontsize = 12,
                     gaps_row = c(4, 10),
                     gaps_col = c(46, 76),
                     width = 8, height = 6, filename = 'figure3_whole_transcriptome_heatmap.pdf')



}

if (sys.nframe() == 0) {
  main()
}
