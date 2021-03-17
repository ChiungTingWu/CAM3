#' Select markers in S matrix with COT
#'
#' This function computes cosine values to select markers.
#' @param data A data set that will be internally coerced into a matrix.
#'     If data is not NULL, "thres.low" and "thres.high" would decide which
#'     genes should be removed in data. If data is NULL, the thresholds work
#'     on "Sest" directly.
#' @param Sest Subpopulation-specific expression matrix. Each row is a gene 
#'     and each column is a sample.
#' @param thres.low The lower bound of percentage of genes to keep with ranked
#'     norm. The value should be between 0 and 1. The default is 0.05.
#' @param thres.high The higher bound of percentage of genes to keep with ranked
#'     norm. The value should be between 0 and 1. The default is 1.
#' @param cos.thres The cosine threshold for markers.
#' @param top The upper bound of total markers.
#' @param per The upper bound of markers of each Subpopulation.
#' @details This function uses COT to detect markers.
#' @return A list of vectors, each of which contains marker genes for one
#' subpopulation.
#' @export
#' @examples
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #Find marker genes
#' MGlist <- cotMG(data, rASest3$S, cos.thres = 0.99)
cotMG <- function(data=NULL, Sest, thres.low=0.05, thres.high=1, cos.thres=1, 
                  top=NULL, per=NULL) {
    
    if (thres.low >= thres.high) {
        stop("thres.low must be smaller than thres.high!")
    }
    Valid <- NULL
    nonzero.ind <- NULL
    if (!is.null(data)) {
        if (is(data, "data.frame")) {
            data <- as.matrix(data)
        } else if (is(data, "SummarizedExperiment")) {
            data <- SummarizedExperiment::assay(data)
        } else if (is(data, "ExpressionSet")) {
            data <- Biobase::exprs(data)
        } else if (is(data, "matrix") == FALSE) {
            stop("Only matrix, dataframe, SummarizedExperiment and ExpressionSet
            object are supported for expression data!")
        }
        if (nrow(data) != nrow(Sest)){
            stop("The row number of data and S matrix should be the same!")
        }
        if (sum(is.na(data)) > 0) {
            stop("Data with missing values are not supported!")
        }
        if (sum(data<0) > 0) {
            stop("Only non-negative data are supported!")
        }
        if (is.null(rownames(data))) {
            rownames(data) <- seq_len(nrow(data))
        }
        nonzero.ind <- rowSums(data) > 0
        data <- data[nonzero.ind, ]
        X <- t(data)
        
        sigNorm <- apply(X, 2, function(x) norm(matrix(x),"F") )
        Valid <- sigNorm >= quantile(sigNorm, thres.low) &
            sigNorm <= quantile(sigNorm, thres.high)
        X <- X[, Valid]
    }
    
    if (is(Sest, "data.frame")) {
        Sest <- as.matrix(Sest)
    } else if (is(Sest, "matrix") == FALSE) {
        stop("Only matrix and data frame are supported for S matrix!")
    }
    if (sum(is.na(Sest)) > 0) {
        stop("S matrix with missing values are not supported!")
    }
    if (sum(Sest<0) > 0) {
        stop("Only non-negative S matrix are supported!")
    }
    if (!is.null(nonzero.ind)){
        Sest <- Sest[nonzero.ind, ]
    }
    if (is.null(rownames(Sest))) {
        rownames(Sest) <- seq_len(nrow(Sest))
    }
    S <- t(Sest)
    if (is.null(Valid)) {
        sigNorm <- apply(S, 2, function(x) norm(matrix(x),"F") )
        Valid <- sigNorm >= quantile(sigNorm, thres.low) &
            sigNorm <= quantile(sigNorm, thres.high)
    }
    S <- S[, Valid]
    S.norm <- apply(S, 2, function(x) x / norm(matrix(x),"2"))
    mg.cos <- matrix(NA, ncol(S), 2)
    rownames(mg.cos) <- colnames(S)
    colnames(mg.cos) <- c('cos', 'source')
    for (i in seq_len(ncol(S))){
        mg.cos[i, 1] <- max(S.norm[, i])
        mg.cos[i, 2] <- which.max(S.norm[, i])
    }
    
    if (!is.null(cos.thres)){
        mg.cos <- mg.cos[mg.cos[, 1] >= cos.thres, ]
    }
    
    mg.cos.sorted <- mg.cos[sort.int(mg.cos[, 1], decreasing=TRUE, 
                                     index.return=TRUE)$ix, ]
    if (!is.null(top)){
        mg.cos.sorted <- mg.cos.sorted[1:top, ]
    }
        
    mg.list <- vector("list", nrow(S))
    for (i in seq_along(mg.list)){
        mg.list[[i]] <- rownames(mg.cos.sorted[mg.cos.sorted[, 2] == i, ])
    }
    
    if (!is.null(per)){
        for (i in seq_along(mg.list)){
            if (length(mg.list[[i]]) > per){
                mg.list[[i]] <- mg.list[[i]][1:per]
            }
        }
    }
    
    return(list(mg.list = mg.list, mg.cos = mg.cos))
}
