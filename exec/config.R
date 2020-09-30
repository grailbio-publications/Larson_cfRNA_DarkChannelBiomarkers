# set packrat library path
.libPaths(packrat::lib_dir('../'))

if (!is.na(Sys.getenv("R_STRICT_WARNINGS", unset = NA))) {
  options(warn = 1)
  options(warnPartialMatchArgs = TRUE)
  options(warnPartialMatchAttr = TRUE)
  options(warnPartialMatchDollar = TRUE)
}


library(cellfreetranscriptome)

`%>%` <- magrittr::`%>%`

pilot_breast_dcbs <- c('CSN1S1', 'FABP7', 'OPN1SW', 'SCGB2A2',
                       'LALBA', 'CASP14', 'KLK5', 'WFDC2')
pilot_lung_dcbs <- c('SLC34A2', 'GABRG1', 'ROS1', 'AGR2', 'GNAT3',
                     'SFTPA2', 'MUC5B', 'SFTA3', 'SMIM22', 'CXCL17',
                     'BPIFA1', 'WFDC2')
heteroDE <- c('SCGB2A2', 'CASP14', 'FABP7', 'CRABP2', 'VGLL1', 'SERPINB5', 'TFF1')
biscotti_pos_control_genes <- c('ALK', 'FGFR4', 'NTRK1', 'EGFR', 'FGFR3', 'FGFR2',
                                'MET', 'RET', 'FGFR1', 'ERBB2', 'KRAS', 'PIK3CA')
