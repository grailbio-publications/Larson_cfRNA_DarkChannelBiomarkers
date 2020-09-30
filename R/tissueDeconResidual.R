#' Calculate tissue specificity score (TSS)
#'
#' Calculates tissue specificity score (TSS) from median tissue expression matrix,
#' the fraction of each tissue expression over the sum of all tissue expression
#'
#' @param tissue_matrix data_frame with columns of gene-wise median counts by reference tissue
#' @param background_threshold numeric value at which to binarize tissue-specific expression
#'
#' @return A data_frame with one row per containing the tissue specificity score for each tissue
#' @export
#'
compute_TSS <- function(tissue_matrix, background_threshold=1) {
  # expr_weight is the tissue expression over the sum of all tissue expression
  tissue_matrix <- tissue_matrix+1 # add 1 to avoid log 0 when calculating shannon enthopy
  expr_weight_matrix <- as.data.frame(matrix(nrow=nrow(tissue_matrix),
                                             ncol=ncol(tissue_matrix)))
  # identify the dominant tissue of each gene matrix
  TSS_matrix_GTEx <- as.data.frame(matrix(nrow=nrow(tissue_matrix),
                                          ncol=7))
  colnames(TSS_matrix_GTEx) <- c("gene_id", "Top_tissue", "Tissue_expr", "Expr_weight",
                                 "Entropy", "NumOfTissue_on", "Tissue_on")

  # calculate the number tissue the gene expresses
  num_of_tissue_with_gene_expression<-rep(NA, nrow(tissue_matrix))
  names(num_of_tissue_with_gene_expression)<-rownames(tissue_matrix)
  for (i in seq_len(nrow(tissue_matrix))) {
    # gene expression vector across all tissue
    gene_tissue_vec <- tissue_matrix[i, ]
    # tissue expression fraction at each tissue
    TSS_weight <- gene_tissue_vec/sum(gene_tissue_vec)

    TSS_matrix_GTEx$gene_id[i] <- rownames(tissue_matrix)[i]
    TSS_matrix_GTEx$Top_tissue[i] <- colnames(tissue_matrix)[order(TSS_weight,
                                                                   decreasing=T)[1]]
    TSS_matrix_GTEx$Tissue_expr[i] <- gene_tissue_vec[order(TSS_weight,
                                                            decreasing=T)[1]]
    TSS_matrix_GTEx$Expr_weight[i] <- as.numeric(TSS_weight[order(TSS_weight,
                                                                  decreasing=T)[1]])
    # calculate Shannon entropy
    TSS_matrix_GTEx$Entropy[i] <- sum(sapply(TSS_weight,
                                             function(x) -x*log2(x)))
    # calculate the number of the tissue this gene expresses
    TSS_matrix_GTEx$NumOfTissue_on[i] <- table(gene_tissue_vec > background_threshold)["TRUE"]
    if (is.na(TSS_matrix_GTEx$NumOfTissue_on[i])) {
      TSS_matrix_GTEx$NumOfTissue_on[i] <- 0
    } else if (TSS_matrix_GTEx$NumOfTissue_on[i]<6) {
      TSS_matrix_GTEx$Tissue_on[i] <- paste(colnames(tissue_matrix)[gene_tissue_vec >
                                                                      background_threshold],
                                            collapse = ",")
    }
    #the fraction of each tissue expression over the sum of all tissue expression
    expr_weight_matrix[i, ] <- TSS_weight
  }
  return(TSS_matrix_GTEx)
}


#' Estimates tissue deconvolution using quadprog
#'
#' @param signatures data_frame with columns containing tissue-specific gene of different tissue.
#' @param dataset_DeconRNASeq data_frame with columns of mixture RNA-seq samples.
#'
#' @return A data_frame with one row per containing the estimation of tissue fraction
#' @export
#'
tissueDeconResidual <- function(signatures, dataset_DeconRNASeq) {
  out.all <- c()
  for (k in seq_len(ncol(dataset_DeconRNASeq))) {
    yvec<- as.vector(dataset_DeconRNASeq[, k])

    # Only require each cell type contribution to be >= 0 and their sum to be
    # <= 1.
    basis <- as.matrix(signatures)
    d <- ncol(basis)
    m <- nrow(basis)
    Amat <- t(rbind(matrix(rep(-1, d), nrow = 1), diag(d)))
    bvec <- c(-1, rep(0, d))

    # Setup cross-product matrix and vector for optimization.
    Dmat <- crossprod(basis / max(basis))
    dvec <- crossprod(basis / max(basis), yvec / max(basis))

    # Solve the quadratic program and return a tidy form of the result.
    solution <- quadprog::solve.QP(Dmat, dvec, Amat, bvec)
    proportions <- solution$solution
    remainder <- 1 - sum(proportions)
    result <- c(proportions, remainder)
    out.all <- rbind(out.all, result)
  }

  fraction_predicted <- round(out.all, digits=6)
  rownames(fraction_predicted) <- colnames(dataset_DeconRNASeq)
  colnames(fraction_predicted) <- c(colnames(signatures), "Remainder")

  return(fraction_predicted)
}
