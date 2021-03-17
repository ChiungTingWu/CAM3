#' A and S matrix estimation by CAM3
#'
#' This function estimates A and S matrix based on marker gene clusters
#' detected by CAM3.
#' @param MGResult An object of class "\code{\link[debCAM]{CAMMGObj}}" obtained 
#'     from \code{\link{CAM3MGCluster}} function.
#' @param PrepResult An object of class "\code{\link[debCAM]{CAMPrepObj}}" 
#'     obtained from \code{\link[debCAM]{CAMPrep}} function.
#' @param data Matrix of mixture expression profiles which need to be the same
#'     as the input of \code{\link[debCAM]{CAMPrep}}.
#'     Data frame, SummarizedExperiment or ExpressionSet object will be
#'     internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     Data should be in non-log linear space with non-negative numerical values
#'     (i.e. >= 0). Missing values are not supported.
#'     All-zero rows will be removed internally.
#' @param corner.strategy The method to detect corner clusters.
#'     1: minimum sum of margin-of-errors; 2: minimum sum of reconstruction
#'     errors. The default is 2.
#' @param generalNMF If TRUE, the decomposed proportion matrix has no sum-to-one
#'     constraint for each row. Without this constraint, the scale ambiguity of
#'     each column vector in proportion matrix will not be removed.
#'     The default is FALSE.
#' @details This function is used internally by \code{\link{CAM3Run}} function to
#' estimate proportion matrix (A), subpopulation-specific expression matrix (S)
#' and mdl values. It can also be used when you want to perform CAM step by
#' step.
#'
#' The mdl value is based on original data with A matrix estimated by 
#' transforming dimension-reduced A matrix back to original space.
#' The mdl value is the sum of two terms: code length of data under the model
#' and code length of model. Both mdl value and the first term (code length
#' of data) will be returned.
#' @return An object of class "\code{\link[debCAM]{CAMASObj}}" containing the
#' following components:
#' \item{Aest}{Estimated proportion matrix.}
#' \item{Sest}{Estimated subpopulation-specific expression matrix.}
#' \item{Aest.proj}{Estimated proportion matrix, before removing scale 
#'     ambiguity.}
#' \item{Ascale}{The estimated scales to remove scale ambiguity
#'     of each column vector in Aest. Sum-to-one constraint on each row of
#'     Aest is used for scale estimation.}
#' \item{AestO}{No use in CAM3, will be deprecated in future version.}
#' \item{SestO}{No use in CAM3, will be deprecated in future version.}
#' \item{AestO.proj}{No use in CAM3, will be deprecated in future version.}
#' \item{AscaleO}{No use in CAM3, will be deprecated in future version.}
#' \item{datalength}{Code length of data.}
#' \item{mdl}{The mdl value.}
#' @export
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #preprocess data
#' rPrep3 <- CAMPrep3(data, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
#'
#' #Marker gene cluster detection with greedy search
#' rMGC3 <- CAM3MGCluster(3, rPrep3)
#'
#' #A and S matrix estimation
#' rASest3 <- CAM3ASest(rMGC3, rPrep3, data)
CAM3ASest <- function(MGResult, PrepResult, data, corner.strategy = 2, 
                    generalNMF = FALSE) {
    if (is.null(MGResult)) {
        return (NULL)
    }
    if (is(data, "data.frame")) {
        data <- as.matrix(data)
    } else if (is(data, "SummarizedExperiment")) {
        data <- SummarizedExperiment::assay(data)
    } else if (is(data, "ExpressionSet")) {
        data <- Biobase::exprs(data)
    } else if (is(data, "matrix") == FALSE) {
        stop("Only matrix, data frame, SummarizedExperiment and ExpressionSet
            object are supported for expression data!")
    }
    if (is.null(rownames(data))) {
        rownames(data) <- seq_len(nrow(data))
        warning('Gene/probe name is missing!')
    }
    data <- data[rowSums(data) > 0,]

    c.valid <- !(PrepResult@cluster$cluster %in% PrepResult@c.outlier)
    geneValid <- PrepResult@Valid
    geneValid[geneValid][!c.valid] <- FALSE

    Xprep <- PrepResult@Xprep[,c.valid]
    Xproj <- PrepResult@Xproj[,c.valid]
    dataSize <- ncol(Xprep)

    Aproj <- PrepResult@centers[,as.character(MGResult@corner[corner.strategy,])]
    Kest <- ncol(Aproj)

    scale <- rep(1, Kest)
    if (generalNMF == FALSE) {
        scale <- tryCatch({
            c(NMF::.fcnnls(Aproj, PrepResult@W%*%PrepResult@SW)$coef)
        }, error = function(msg) {
            res <- apply(PrepResult@W%*%PrepResult@SW, 2, function(x) (nnls::nnls(Aproj,x))$x)
            return(as.vector(res))
        })
        scale[scale<1e-10] <- 0.01/(sqrt(colSums(Aproj^2)))[scale<1e-10]
    }
    Apca <- Aproj %*% diag(scale)


    if (ncol(PrepResult@W) != ncol(data)) {
        stop("Data should be the same with the input of CAMPrep()!")
    }
    X <- t(data) * PrepResult@SW

    Aproj.org <- corpcor::pseudoinverse(PrepResult@W) %*% Aproj
    scale.org <- rep(1, Kest)
    if (generalNMF == FALSE) {
        scale.org <- tryCatch({
            c(NMF::.fcnnls(Aproj.org, matrix(PrepResult@SW, ncol=1))$coef)
        }, error = function(msg) {
            res <- apply(matrix(PrepResult@SW, ncol=1), 2, function(x) (nnls::nnls(Aproj.org,x))$x)
            return(as.vector(res))
        })
        scale.org[scale.org<1e-10] <-
            0.01/(sqrt(colSums(Aproj.org^2)))[scale.org<1e-10]
    }
    Aest.org <- Aproj.org %*% diag(scale.org)
    Aest.org[Aest.org < 0] <- 0
    S2 <- tryCatch({
        NMF::.fcnnls(Aest.org, X)$coef
    }, error = function(msg) {
        return(apply(X, 2, function(x) (nnls::nnls(Aest.org,x))$x))
    })
    datalength2 <- (nrow(X) * dataSize) / 2 *
        log(mean((as.vector(X[,geneValid] - Aest.org %*% S2[,geneValid]))^2))
    reconErr2 <- sum((as.vector(X[,geneValid] - Aest.org %*% S2[,geneValid]))^2)
    penalty2 <- ((Kest - 1) * nrow(X)) / 2 * log(dataSize) +
        (Kest * dataSize) / 2 * log(nrow(X))
    mdl2 <- datalength2 + penalty2
    if (generalNMF == FALSE) {
        Aest.org <- Aest.org/rowSums(Aest.org)
    }

    return(new("CAMASObj", 
                Aest=Aest.org, Sest=t(S2), Aest.proj=Aproj.org, Ascale=scale.org,
                AestO=Aest.org, SestO=t(S2), AestO.proj=Aproj.org, 
                AscaleO=scale.org, datalength=datalength2, mdl=mdl2))
}
