source('./config.R')

parse_arguments <- function() {

  parser <- argparse::ArgumentParser()

  parser$add_argument('--demographics',
                      default = 'sample_demographics.tsv',
                      required = TRUE,
                      help = 'TSV file with CCGA sample cohort demographics')
  parser$add_argument('--count_cfrna',
                      required = TRUE,
                      default = 'whole_transcriptome_cfrna_spliced_counts.tsv',
                      help = 'TSV file with whole transcriptome plasma gene counts')
  parser$add_argument('--count_tissue',
                      default = 'whole_transcriptome_tissue_spliced_counts.tsv',
                      required = TRUE,
                      help = '')
  parser$add_argument('--tcga_luad',
                      required = TRUE,
                      default = 'luad_tumor_counts.tsv',
                      help = '')
  parser$add_argument('--tcga_lusc',
                      required = TRUE,
                      default = 'lusc_tumor_counts.tsv',
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

  sample_metadata <- readr::read_tsv(args$demographics) %>%
    dplyr::mutate(primccat = dplyr::case_when(grepl('BrCa', `Sample ID`) ~ 'Breast',
                                              grepl('LuCa', `Sample ID`) ~ 'Lung',
                                              grepl('NC', `Sample ID`) ~ 'Non-cancer'),
                  Subtype = dplyr::case_when(primccat == 'Non-cancer' ~ 'Non-cancer',
                                             TRUE ~ Subtype))
  # get ccga plasma counts
  count_cfrna_df <- readr::read_tsv(args$count_cfrna)
  # get ccga tissue counts
  count_tissue_df <- readr::read_tsv(args$count_tissue)

  # get tcga tissue counts
  tcga_luad_tumor_counts <- readr::read_tsv(args$tcga_luad)
  tcga_lusc_tumor_counts <- readr::read_tsv(args$tcga_lusc)
  tcga_hr_pos_her2_neg_counts <- readr::read_tsv(args$tcga_hrpos_her2neg)
  tcga_hr_pos_her2_pos_counts <- readr::read_tsv(args$tcga_hrpos_her2pos)
  tcga_hr_neg_her2_pos_counts <- readr::read_tsv(args$tcga_hrneg_her2pos)
  tcga_tnbc_counts <- readr::read_tsv(args$tcga_tnbc)

  # sample RPMs
  rpm_plasma <- get_rpm(count_cfrna_df)
  rpm_tissue <- get_rpm(count_tissue_df)
  tcga_hr_pos_her2_neg_rpm <- get_rpm(tcga_hr_pos_her2_neg_counts) %>%
    tidyr::gather(`Sample ID`, 'TCGA expr', -gene_id)
  tcga_hr_pos_her2_pos_rpm <- get_rpm(tcga_hr_pos_her2_pos_counts) %>%
    tidyr::gather(`Sample ID`, 'TCGA expr', -gene_id)
  tcga_hr_neg_her2_pos_rpm <- get_rpm(tcga_hr_neg_her2_pos_counts) %>%
    tidyr::gather(`Sample ID`, 'TCGA expr', -gene_id)
  tcga_tnbc_rpm <- get_rpm(tcga_tnbc_counts) %>%
    tidyr::gather(`Sample ID`, 'TCGA expr', -gene_id)
  tcga_luad_tumor_rpm <- get_rpm(tcga_luad_tumor_counts) %>%
    tidyr::gather(`Sample ID`, 'TCGA expr', -gene_id)
  tcga_lusc_tumor_rpm <- get_rpm(tcga_lusc_tumor_counts) %>%
    tidyr::gather(`Sample ID`, 'TCGA expr', -gene_id)

  #----------------------------------------------------------------------#
  #----------------CCGA breast cancer plasma plot------------------------#
  #----------------------------------------------------------------------#
  # breast cancer subtypes are annotated as
  # TNBC (HR-/HER2-), HR+ (HR+/HER2-, HR+/HER2), or HR- (HR-/HER2+)
  d2plot <- count_cfrna_df %>%
    tidyr::gather(`Sample ID`, 'cfRNA expr', -gene_id) %>%
    dplyr::left_join(sample_metadata) %>%
    dplyr::filter(gene_id %in% c('FABP7', 'SCGB2A2')) %>%
    dplyr::mutate(Subtype = dplyr::case_when(Subtype == "HR+/HER2+" ~ 'HR+',
                                             Subtype == "HR+/HER2-" ~ 'HR+',
                                             Subtype == "HR-/HER2+" ~ 'HR-',
                                             TRUE ~ Subtype)) %>%
    dplyr::filter(Subtype %in% c('Non-cancer',
                                 'HR+',
                                 'TNBC')) %>%
    dplyr::mutate(Subtype = factor(Subtype,
                                   levels = c('Non-cancer', 'HR+', 'TNBC')),
                  expression_new = `cfRNA expr`)

  # wilcoxon test - test for differential expression between subtypes
  test_fabp7 <- data.frame(group1 = c('TNBC', 'TNBC'),
                           group2 = c('Non-cancer', 'HR+'),
                           gene = 'FABP7')
  test_fabp7$p.value <- apply(test_fabp7,
                              1,
                              function(x) wilcoxon_test(d2plot,
                                                        x['gene'],
                                                        unlist(x['group1']),
                                                        unlist(x['group2'])))
  test_fabp7 <- test_fabp7 %>%
    dplyr::mutate(gene_id = 'FABP7',
                  sig.level = set_sig_level(p.value),
                  y.position = c(1.4, 1.2))

  test_scgb2a2 <- data.frame(group1 = c('HR+', 'HR+'),
                             group2 = c('Non-cancer', 'TNBC'),
                             gene = 'SCGB2A2')
  test_scgb2a2$p.value <- apply(test_scgb2a2,
                                1,
                                function(x) wilcoxon_test(d2plot,
                                                          x['gene'],
                                                          unlist(x['group1']),
                                                          unlist(x['group2'])))
  test_scgb2a2 <- test_scgb2a2 %>%
    dplyr::mutate(gene_id = 'SCGB2A2',
                  sig.level = set_sig_level(p.value),
                  y.position = c(1.3, 1.5))

  test_dat <- rbind(test_scgb2a2, test_fabp7)

  dummy <- data.frame(gene_id = c('FABP7', 'SCGB2A2', 'FABP7', 'SCGB2A2'),
                      `cfRNA expr` = c(7e-2, 7e-2, 40, 40),
                      Subtype = 'TNBC',
                      stringsAsFactors=FALSE, check.names = F)

  p1_counts_log <- d2plot %>%
    # transforms counts by + 0.5
    dplyr::mutate(`cfRNA expr` = `cfRNA expr` + 0.1) %>%
    ggplot2::ggplot(ggplot2::aes(x = Subtype, y = `cfRNA expr`,
                                 colour = Subtype)) +
    ggplot2::scale_colour_manual(values = c("grey50", "orange1", "navyblue", "red2")) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.25, height = 0) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~gene_id, scales = 'free_y') +
    ggplot2::labs(x = '', y = 'cfRNA expr (count)') +
    ggplot2::scale_y_log10() +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggpubr::stat_pvalue_manual(test_dat, label = "sig.level", linetype = 'dashed',
                               bracket.size = 0.5, label.size = 10) +
    # geom_blank adds padding around axis limits
    ggplot2::geom_blank(data=dummy)

  #----------------------------------------------------------------------#
  #----------------CCGA breast cancer tissue plot------------------------#
  #----------------------------------------------------------------------#
  d2plot <- rpm_tissue %>%
    tidyr::gather(`Sample ID`, 'Tissue expr', -gene_id) %>%
    dplyr::left_join(sample_metadata) %>%
    dplyr::filter(gene_id %in% c('FABP7', 'SCGB2A2')) %>%
    dplyr::mutate(Subtype = dplyr::case_when(Subtype == "HR+/HER2+" ~ 'HR+',
                                             Subtype == "HR+/HER2-" ~ 'HR+',
                                             Subtype == "HR-/HER2+" ~ 'HR-',
                                             TRUE ~ Subtype)) %>%
    dplyr::filter(Subtype %in% c('HR+', 'TNBC')) %>%
    dplyr::mutate(Subtype = factor(Subtype,
                                   levels = c('Non-cancer', 'HR+',
                                              'TNBC')),
                  expression_new = `Tissue expr`)

  # wilcoxon test - test for differential expression between subtypes
  test_fabp7 <- data.frame(group1 = c('TNBC'),
                           group2 = c('HR+'),
                           gene = 'FABP7')
  test_fabp7$p.value <- apply(test_fabp7,
                              1,
                              function(x) wilcoxon_test(d2plot,
                                                        x['gene'],
                                                        unlist(x['group1']),
                                                        unlist(x['group2'])))
  test_fabp7 <- test_fabp7 %>%
    dplyr::mutate(gene_id = 'FABP7',
                  sig.level = set_sig_level(p.value),
                  y.position = c(3.5))

  test_scgb2a2 <- data.frame(group1 = c('HR+'),
                             group2 = c('TNBC'),
                             gene = 'SCGB2A2')
  test_scgb2a2$p.value <- apply(test_scgb2a2,
                                1,
                                function(x) wilcoxon_test(d2plot,
                                                          x['gene'],
                                                          unlist(x['group1']),
                                                          unlist(x['group2'])))
  test_scgb2a2 <- test_scgb2a2 %>%
    dplyr::mutate(gene_id = 'SCGB2A2',
                  sig.level = set_sig_level(p.value),
                  y.position = 4)

  test_dat <- rbind(test_scgb2a2, test_fabp7)

  p2_rpms_log <- d2plot %>%
    ggplot2::ggplot(ggplot2::aes(x = Subtype, y = `Tissue expr`,
                                 colour = Subtype)) +
    ggplot2::scale_colour_manual(values = c("orange1", "navyblue", "red2")) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.25, height = 0) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~gene_id, scales = 'free_y') +
    ggplot2::labs(x = '', y = 'Tissue expr (RPM)') +
    ggplot2::scale_y_log10() +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::scale_x_discrete(limits = c('Non-cancer', 'HR+', 'TNBC')) +
    ggpubr::stat_pvalue_manual(test_dat, label = "sig.level", linetype = 'dashed',
                               bracket.size = 0.5, label.size = 10)

  #----------------------------------------------------------------------#
  #----------------CCGA lung cancer plasma plot------------------------#
  #----------------------------------------------------------------------#
  d2plot <- count_cfrna_df %>%
    tidyr::gather(`Sample ID`, 'cfRNA expr', -gene_id) %>%
    dplyr::left_join(sample_metadata) %>%
    dplyr::filter(gene_id %in% c('SLC34A2', 'SFTPA2', 'CXCL17', 'SFTA3'),
                  primccat %in% c('Non-cancer', 'Lung'),
                  Subtype %in% c('Non-cancer', 'Adenocarcinoma', 'Squamous cell carcinoma')) %>%
    dplyr::mutate(Subtype = factor(Subtype,
                                   levels = c('Non-cancer',
                                              'Squamous cell carcinoma',
                                              'Adenocarcinoma')),
                  subtype_new = Subtype,
                  expression_new = `cfRNA expr`)

  # wilcoxon test - test for differential expression between subtypes
  test_slc34a2 <- data.frame(group1 = c('Adenocarcinoma', 'Adenocarcinoma'),
                             group2 = c('Non-cancer', 'Squamous cell carcinoma'),
                             gene = 'SLC34A2')
  test_slc34a2$p.value <- apply(test_slc34a2,
                                1,
                                function(x) wilcoxon_test(d2plot,
                                                          x['gene'],
                                                          unlist(x['group1']),
                                                          unlist(x['group2'])))
  test_slc34a2 <- test_slc34a2 %>%
    dplyr::mutate(gene_id = 'SLC34A2',
                  sig.level = set_sig_level(p.value),
                  y.position = c(3, 2.7))

  test_sfta3 <- data.frame(group1 = c('Adenocarcinoma', 'Adenocarcinoma'),
                           group2 = c('Non-cancer', 'Squamous cell carcinoma'),
                           gene = 'SFTA3')
  test_sfta3$p.value <- apply(test_sfta3,
                              1,
                              function(x) wilcoxon_test(d2plot,
                                                        x['gene'],
                                                        unlist(x['group1']),
                                                        unlist(x['group2'])))
  test_sfta3 <- test_sfta3 %>%
    dplyr::mutate(gene_id = 'SFTA3',
                  sig.level = set_sig_level(p.value),
                  y.position = c(1.8, 1.6))

  test_cxcl17 <- data.frame(group1 = c('Adenocarcinoma', 'Adenocarcinoma'),
                            group2 = c('Non-cancer', 'Squamous cell carcinoma'),
                            gene = 'CXCL17')
  test_cxcl17$p.value <- apply(test_cxcl17,
                               1,
                               function(x) wilcoxon_test(d2plot,
                                                         x['gene'],
                                                         unlist(x['group1']),
                                                         unlist(x['group2'])))
  test_cxcl17 <- test_cxcl17 %>%
    dplyr::mutate(gene_id = 'CXCL17',
                  sig.level = set_sig_level(p.value),
                  y.position = c(1.60, 1.50))

  test_sftpa2 <- data.frame(group1 = c('Adenocarcinoma', 'Adenocarcinoma'),
                            group2 = c('Non-cancer', 'Squamous cell carcinoma'),
                            gene = 'SFTPA2')
  test_sftpa2$p.value <- apply(test_sftpa2,
                               1,
                               function(x) wilcoxon_test(d2plot,
                                                         x['gene'],
                                                         unlist(x['group1']),
                                                         unlist(x['group2'])))
  test_sftpa2 <- test_sftpa2 %>%
    dplyr::mutate(gene_id = 'SFTPA2',
                  sig.level = set_sig_level(p.value),
                  y.position = c(2.2, 2))

  test_dat <- rbind(test_cxcl17, test_sftpa2) %>%
    rbind(., test_slc34a2) %>%
    rbind(., test_sfta3)

  dummy <- data.frame(gene_id = c('CXCL17', 'SFTPA2', 'SFTA3', 'SLC34A2',
                                  'CXCL17', 'SFTPA2', 'SFTA3', 'SLC34A2'),
                      `cfRNA expr` = c(7e-2, 7e-2, 7e-2, 7e-2,
                                       80, 1e2, 90, 1e3),
                      Subtype = 'Adenocarcinoma',
                      stringsAsFactors=FALSE, check.names = F)

  p3_counts_log <- d2plot %>%
    # transforms counts
    dplyr::mutate(`cfRNA expr` = `cfRNA expr` + 0.1) %>%
    ggplot2::ggplot(ggplot2::aes(x = Subtype, y = `cfRNA expr`, colour = Subtype)) +
    ggplot2::scale_colour_manual(values = c("grey50", "cyan3", "magenta4")) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.25, height = 0) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~gene_id, nrow = 1, scales = 'free_y') +
    ggplot2::labs(x = '', y = 'cfRNA expr (count)') +
    ggplot2::scale_y_log10() +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggpubr::stat_pvalue_manual(test_dat, label = "sig.level", linetype = 'dashed',
                               bracket.size = 0.5, label.size = 10) +
    # geom_blank adds padding around axis limits
    ggplot2::geom_blank(data=dummy)

  #--------------------------------------------------------------------#
  #----------------CCGA lung cancer tissue plot------------------------#
  #--------------------------------------------------------------------#
  d2plot <- rpm_tissue %>%
    tidyr::gather(`Sample ID`, 'Tissue expr', -gene_id) %>%
    dplyr::left_join(sample_metadata) %>%
    dplyr::filter(gene_id %in% c('SLC34A2', 'SFTPA2', 'CXCL17', 'SFTA3'),
                  primccat %in% c('Non-cancer', 'Lung'),
                  Subtype %in% c('Non-cancer', 'Squamous cell carcinoma', 'Adenocarcinoma')) %>%
    dplyr::mutate(Subtype = factor(Subtype,
                                   levels = c('Non-cancer', 'Squamous cell carcinoma',
                                              'Adenocarcinoma')),
                  subtype_new = Subtype,
                  expression_new = `Tissue expr`)

  # wilcoxon test - test for differential expression between subtypes
  test_slc34a2 <- data.frame(group1 = c('Adenocarcinoma'),
                             group2 = c('Squamous cell carcinoma'),
                             gene = 'SLC34A2')
  test_slc34a2$p.value <- apply(test_slc34a2,
                                1,
                                function(x) wilcoxon_test(d2plot,
                                                          x['gene'],
                                                          unlist(x['group1']),
                                                          unlist(x['group2'])))
  test_slc34a2 <- test_slc34a2 %>%
    dplyr::mutate(gene_id = 'SLC34A2',
                  sig.level = set_sig_level(p.value),
                  y.position = 1)

  test_sfta3 <- data.frame(group1 = c('Adenocarcinoma'),
                           group2 = c('Squamous cell carcinoma'),
                           gene = 'SFTA3')
  test_sfta3$p.value <- apply(test_sfta3,
                              1,
                              function(x) wilcoxon_test(d2plot,
                                                        x['gene'],
                                                        unlist(x['group1']),
                                                        unlist(x['group2'])))
  test_sfta3 <- test_sfta3 %>%
    dplyr::mutate(gene_id = 'SFTA3',
                  sig.level = set_sig_level(p.value),
                  y.position = 1)

  test_cxcl17 <- data.frame(group1 = c('Adenocarcinoma'),
                            group2 = c('Squamous cell carcinoma'),
                            gene = 'CXCL17')
  test_cxcl17$p.value <- apply(test_cxcl17,
                               1,
                               function(x) wilcoxon_test(d2plot,
                                                         x['gene'],
                                                         unlist(x['group1']),
                                                         unlist(x['group2'])))
  test_cxcl17 <- test_cxcl17 %>%
    dplyr::mutate(gene_id = 'CXCL17',
                  sig.level = set_sig_level(p.value),
                  y.position = 1)

  test_sftpa2 <- data.frame(group1 = c('Adenocarcinoma'),
                            group2 = c('Squamous cell carcinoma'),
                            gene = 'SFTPA2')
  test_sftpa2$p.value <- apply(test_sftpa2,
                               1,
                               function(x) wilcoxon_test(d2plot,
                                                         x['gene'],
                                                         unlist(x['group1']),
                                                         unlist(x['group2'])))
  test_sftpa2 <- test_sftpa2 %>%
    dplyr::mutate(gene_id = 'SFTPA2',
                  sig.level = set_sig_level(p.value),
                  y.position = 1)

  test_dat <- rbind(test_cxcl17, test_sftpa2) %>%
    rbind(., test_slc34a2) %>%
    rbind(., test_sfta3)

  p4_rpms_log <- d2plot %>%
    ggplot2::ggplot(ggplot2::aes(x = Subtype, y = `Tissue expr`, colour = Subtype)) +
    ggplot2::scale_colour_manual(values = c("cyan3", "magenta4")) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.25, height = 0) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~gene_id, nrow = 1, scales = 'free_y') +
    ggplot2::labs(x = '', y = 'Tissue expr (RPM)') +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_discrete(limits = c('Non-cancer', 'Squamous cell carcinoma',
                                         'Adenocarcinoma')) +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggpubr::stat_pvalue_manual(test_dat, label = "sig.level", linetype = 'dashed',
                               bracket.size = 0.5, label.size = 10)

  #----------------------------------------------------------------------#
  #----------------TCGA breast cancer tissue plot------------------------#
  #----------------------------------------------------------------------#
  d2plot <- rbind(tcga_hr_pos_her2_neg_rpm,
                  tcga_hr_pos_her2_pos_rpm) %>%
    dplyr::mutate(subtype = 'HR+') %>%
    rbind(tcga_tnbc_rpm %>%
            dplyr::mutate(subtype = 'TNBC')) %>%
    dplyr::filter(gene_id %in% c('FABP7', 'SCGB2A2')) %>%
    dplyr::mutate(subtype = factor(subtype, levels = c('Non-cancer', 'HR+', 'TNBC')),
                  Subtype = subtype,
                  expression_new = `TCGA expr`)

  # wilcoxon test - test for differential expression between subtypes
  test_fabp7 <- data.frame(group1 = c('HR+'),
                           group2 = c('TNBC'),
                           gene = 'FABP7')
  test_fabp7$p.value <- apply(test_fabp7,
                              1,
                              function(x) wilcoxon_test(d2plot,
                                                        x['gene'],
                                                        unlist(x['group1']),
                                                        unlist(x['group2'])))
  test_fabp7 <- test_fabp7 %>%
    dplyr::mutate(gene_id = 'FABP7',
                  sig.level = set_sig_level(p.value),
                  y.position = c(3.5))

  test_scgb2a2 <- data.frame(group1 = c('HR+'),
                             group2 = c('TNBC'),
                             gene = 'SCGB2A2')
  test_scgb2a2$p.value <- apply(test_scgb2a2,
                                1,
                                function(x) wilcoxon_test(d2plot,
                                                          x['gene'],
                                                          unlist(x['group1']),
                                                          unlist(x['group2'])))
  test_scgb2a2 <- test_scgb2a2 %>%
    dplyr::mutate(gene_id = 'SCGB2A2',
                  sig.level = set_sig_level(p.value),
                  y.position = c(4.8))

  test_dat <- rbind(test_scgb2a2, test_fabp7)

  p5_rpms_log <- d2plot %>%
    ggplot2::ggplot(ggplot2::aes(x = subtype, y = `TCGA expr`, colour = subtype)) +
    ggplot2::scale_colour_manual(values = c("orange1", "navyblue", "red2")) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.25, height = 0, alpha = 0.25) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~gene_id, scales = 'free_y') +
    ggplot2::labs(x = '', y = 'TCGA expr (RPM)') +
    ggplot2::scale_y_log10() +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::scale_x_discrete(limits = c('Non-cancer', 'HR+', 'TNBC')) +
    ggpubr::stat_pvalue_manual(test_dat, label = "sig.level", linetype = 'dashed',
                               bracket.size = 0.5, label.size = 10) +
    ggplot2::coord_cartesian(clip = "off")

  #-------------------------------------------------------------------------------------#
  #----------------TCGA breast cancer tissue plot (supplemental)------------------------#
  #-------------------------------------------------------------------------------------#
  d2plot <- rbind(tcga_hr_pos_her2_neg_rpm,
                  tcga_hr_pos_her2_pos_rpm) %>%
    dplyr::mutate(subtype = 'HR+') %>%
    rbind(tcga_tnbc_rpm %>%
            dplyr::mutate(subtype = 'TNBC')) %>%
    dplyr::filter(gene_id %in% pilot_breast_dcbs) %>%
    dplyr::mutate(subtype = factor(subtype, levels = c('Non-cancer', 'HR+', 'TNBC')),
                  Subtype = subtype,
                  expression_new = `TCGA expr`)

  p5_rpms_log_supp <- d2plot %>%
    ggplot2::ggplot(ggplot2::aes(x = subtype, y = `TCGA expr`, colour = subtype)) +
    ggplot2::scale_colour_manual(values = c("orange1", "navyblue", "red2")) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.25, height = 0, alpha = 0.25) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~gene_id, scales = 'free_y', nrow = 2) +
    ggplot2::labs(x = '', y = 'TCGA expr (RPM)') +
    ggplot2::scale_y_log10() +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::scale_x_discrete(limits = c('Non-cancer', 'HR+', 'TNBC')) +
  ggplot2::coord_cartesian(clip = "off")
  p5_rpms_log_supp
  ggplot2::ggsave('figure4_tcga_breast_expression.pdf', device = 'pdf',
                  width = 16, height = 8, dpi = 300)

  #--------------------------------------------------------------------#
  #----------------TCGA lung cancer tissue plot------------------------#
  #--------------------------------------------------------------------#
  d2plot <- rbind(tcga_luad_tumor_rpm %>% dplyr::mutate(subtype = 'Adenocarcinoma'),
                  tcga_lusc_tumor_rpm %>% dplyr::mutate(subtype = 'Squamous cell carcinoma')) %>%
    dplyr::filter(gene_id %in% c('SLC34A2', 'SFTPA2', 'CXCL17', 'SFTA3')) %>%
    dplyr::mutate(subtype = factor(subtype,
                                   levels = c('Non-cancer', 'Squamous cell carcinoma',
                                              'Adenocarcinoma')),
                  Subtype = subtype,
                  expression_new = `TCGA expr`)

  # wilcoxon test - test for differential expression between subtypes
  test_cxcl17 <- wilcox.test(`TCGA expr` ~ subtype, data=d2plot %>%
                               dplyr::filter(gene_id == 'CXCL17'))$p.value %>%
    as.data.frame() %>%
    setNames('p.value') %>%
    dplyr::mutate(gene_id = 'CXCL17',
                  sig.level = set_sig_level(p.value),
                  group1 = 'Adenocarcinoma',
                  group2 = 'Squamous cell carcinoma',
                  y.position = 1e03)
  test_sfta3 <- wilcox.test(`TCGA expr` ~ subtype, data=d2plot %>%
                              dplyr::filter(gene_id == 'SFTA3'))$p.value %>%
    as.data.frame() %>%
    setNames('p.value') %>%
    dplyr::mutate(gene_id = 'SFTA3',
                  sig.level = set_sig_level(p.value),
                  group1 = 'Adenocarcinoma',
                  group2 = 'Squamous cell carcinoma',
                  y.position = 1e03)
  test_sftpa2 <- wilcox.test(`TCGA expr` ~ subtype, data=d2plot %>%
                               dplyr::filter(gene_id == 'SFTPA2'))$p.value %>%
    as.data.frame() %>%
    setNames('p.value') %>%
    dplyr::mutate(gene_id = 'SFTPA2',
                  sig.level = set_sig_level(p.value),
                  group1 = 'Adenocarcinoma',
                  group2 = 'Squamous cell carcinoma',
                  y.position = 1e03)
  test_slc34a2 <- wilcox.test(`TCGA expr` ~ subtype, data=d2plot %>%
                                dplyr::filter(gene_id == 'SLC34A2'))$p.value %>%
    as.data.frame() %>%
    setNames('p.value') %>%
    dplyr::mutate(gene_id = 'SLC34A2',
                  sig.level = set_sig_level(p.value),
                  group1 = 'Adenocarcinoma',
                  group2 = 'Squamous cell carcinoma',
                  y.position = 1e03)

  test_dat <- rbind(test_cxcl17, test_sftpa2) %>%
    rbind(., test_slc34a2) %>%
    rbind(., test_sfta3)

  p6_rpms_log <- d2plot %>%
    #dplyr::mutate(`TCGA expr` = `TCGA expr` + 1e-2) %>%
    ggplot2::ggplot(ggplot2::aes(x = subtype, y = `TCGA expr`, colour = subtype)) +
    ggplot2::scale_colour_manual(values = c("cyan3", "magenta4")) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.25, height = 0, alpha = 0.25) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~gene_id, nrow = 1, scales = 'free_y') +
    ggplot2::labs(x = '', y = 'TCGA expr (RPM)') +
    ggplot2::scale_y_log10() +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::scale_x_discrete(limits = c('Non-cancer', 'Squamous cell carcinoma',
                                         'Adenocarcinoma')) +
    ggpubr::stat_pvalue_manual(test_dat, y.position=5,
                               label = "sig.level", linetype = 'dashed',
                               bracket.size = 0.5, label.size = 10)

  #----------------------------------------------------------------------#
  #----------------TCGA lung cancer tissue plot (supplemental)------------------------#
  #----------------------------------------------------------------------#
  d2plot <- rbind(tcga_luad_tumor_rpm %>% dplyr::mutate(subtype = 'Adenocarcinoma'),
                  tcga_lusc_tumor_rpm %>% dplyr::mutate(subtype = 'Squamous cell carcinoma')) %>%
    dplyr::filter(gene_id %in% pilot_lung_dcbs) %>%
    dplyr::mutate(subtype = factor(subtype,
                                   levels = c('Non-cancer', 'Squamous cell carcinoma',
                                              'Adenocarcinoma')),
                  Subtype = subtype,
                  expression_new = `TCGA expr`)

  p6_rpms_log_supp <- d2plot %>%
    ggplot2::ggplot(ggplot2::aes(x = subtype, y = `TCGA expr`, colour = subtype)) +
    ggplot2::scale_colour_manual(values = c("cyan3", "magenta4")) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.25, height = 0, alpha = 0.25) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~gene_id, scales = 'free_y') +
    ggplot2::labs(x = '', y = 'TCGA expr (RPM)') +
    ggplot2::scale_y_log10() +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::scale_x_discrete(limits = c('Non-cancer', 'Squamous cell carcinoma',
                                         'Adenocarcinoma'))
  p6_rpms_log_supp
  ggplot2::ggsave('figure4_tcga_lung_expression.pdf', device = 'pdf',
                  width = 16, height = 12, dpi = 300)

  #----------------------------------------------------------------------#
  #---------------Arrange figures with cowplot------------------------#
  #----------------------------------------------------------------------#
  p_figs <- cowplot::plot_grid(p1_counts_log +
                                 ggplot2::theme(legend.position="none",
                                                axis.ticks.x = ggplot2::element_blank()),
                               p3_counts_log +
                                 ggplot2::theme(legend.position="none",
                                                axis.ticks.x = ggplot2::element_blank()) +
                                 ggplot2::labs(y=''),
                               p2_rpms_log +
                                 ggplot2::theme(legend.position="none",
                                                axis.ticks.x = ggplot2::element_blank()),
                               p4_rpms_log +
                                 ggplot2::theme(legend.position="none",
                                                axis.ticks.x = ggplot2::element_blank()) +
                                 ggplot2::labs(y=''),
                               p5_rpms_log +
                                 ggplot2::theme(legend.position="none",
                                                axis.ticks.x = ggplot2::element_blank()),
                               p6_rpms_log +
                                 ggplot2::theme(legend.position="none",
                                                axis.ticks.x = ggplot2::element_blank()) +
                                 ggplot2::labs(y=''),
                               ncol=2, rel_widths = c(2, 5),
                               align = 'v', axis = 'l', greedy = F,
                               labels = c('A', 'D', 'B',
                                          'E', 'C', 'F'),
                               label_size = 26)
  p_legend <- cowplot::plot_grid(cowplot::get_legend(
    p1_counts_log +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
      ggplot2::theme(legend.position = "bottom")),
    cowplot::get_legend(
      p3_counts_log +
        ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
        ggplot2::theme(legend.position = "bottom")),
    ncol=2, rel_widths = c(2, 5), align = 'hv')
  p_final <- cowplot::plot_grid(p_figs,
                                p_legend,
                                ncol = 1,
                                rel_heights = c(3, 0.5))
  p_final
  ggplot2::ggsave('figure4_subtype_expression_final.pdf', device = 'pdf',
                  width = 20, height = 12, dpi = 300)



}

if (sys.nframe() == 0) {
  main()
}
