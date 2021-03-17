#' The plot of scatter simplex
#'
#' This function shows scatter simplex of mixture expressions.
#' @param data A data set that will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     Data should be in non-log linear space with non-negative numerical values
#'     (i.e. >= 0). Missing values are not supported.
#'     All-zero rows will be removed internally.
#' @param A Prior/Estimated proportion matrix.
#' @param MGlist A list of vectors, each of which contains known markers
#'     and/or CAM-detected markers for one subpopulation.
#' @param data.extra Extra data points to be shown in the simplex plot, e.g.
#'     the points associated with prior/estimated proportion vectors. The format
#'     should be consistent with \code{data}, so prior/estimated proportion
#'     matrix needs to be transposed before as the input.
#' @param corner.order The order to show simplex corners counterclockwise.
#' @param col The color for data points. The default is "gray".
#' @param pch The symbol/character for data points. The default is 1.
#' @param cex The symbol/character expansion for data points.
#'     The default is 0.8.
#' @param mg.col The colors for marker genes. Marker genes of one subpopulation
#'     could have their own color if a vector is provided. The default is "red".
#' @param mg.pch The symbols/characters for marker genes. Marker genes of one
#'     subpopulation could have their own symbol/character if a vector is
#'     provided. The default is 1.
#' @param mg.cex The symbol/character expansion for marker genes.
#'     The default is 1.2.
#' @param ex.col The colors for extra data points. Each data point could have
#'     its own color if a vector is provided. The default is "black".
#' @param ex.pch The symbols/characters for extra data points. Each data point
#'     could have its own symbol/character if a vector is provided.
#'     The default is 19.
#' @param ex.cex The symbol/character expansion for extra data points.
#'     The default is 1.5.
#' @param ... All other arguments are passed to the plotting command.
#' @details This function can show the scatter simplex and detected marker genes
#' in a 2D plot. The corners in the high-dimensional simplex will still locate
#' at extreme points of low-dimensional simplex.
#' @return A plot to the current device.
#' @export
#' @examples
#' #obtain data, A matrix, marker genes
#' data(ratMix3)
#' data <- ratMix3$X
#' A <- ratMix3$A
#' 
#' #Find marker genes
#' MGlist <- cotMG(data, rASest3$S, cos.thres = 0.99)
#' 
#' #plot simplex for data
#' simplexplot2(data, A)
#' simplexplot2(data, A, MGlist = MGlist) #Color marker genes in the plot
#' simplexplot2(data, A, MGlist = MGlist,
#'     data.extra = t(A)) #show the location of proportion vectors in the plot
#'
#' #set different colors
#' simplexplot2(data, A, MGlist = MGlist,
#'             col = "blue", mg.col = c("red","orange","green"))
simplexplot2 <- function(data, A, MGlist = NULL, data.extra = NULL,
                        col = 'gray', pch = 1, cex = 0.8,
                        mg.col = 'red', mg.pch = 1, mg.cex = 1.2,
                        ex.col = "black", ex.pch = 19, ex.cex = 1.5, ...){
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
    if (sum(data<0) > 0) {
        stop("Only non-negative data are supported!")
    }
    if (nrow(A) != ncol(data)) {
        stop("The number of samples in data and in A matrix
            should be the same!")
    }

    data <- data[rowSums(data) > 0,]
    K <- ncol(A)
    if (is.null(corner.order)) {
        corner.order <- seq_len(K)
    } else if (sum(!(corner.order %in% seq_len(K))) > 0) {
        stop("corner.order must be a permutation of 1 ~", K)
    }

    if(length(mg.col) == 1) {
        mg.col <- rep(mg.col, length(MGlist))
    } else if (length(mg.col) != length(MGlist)) {
        stop("mg.col needs", length(MGlist), "colors!")
    }
    if(length(mg.pch) == 1) {
        mg.pch <- rep(mg.pch, length(MGlist))
    } else if (length(mg.pch) != length(MGlist)) {
        stop("mg.pch needs", length(MGlist), "values!")
    }

    if (!is.null(data.extra)) {
        if (!(is(data.extra, "matrix") | is(data.extra, "data.frame"))) {
            stop("data.extra must be a matrix or a data frame!")
        }
        if (nrow(A) != ncol(data.extra)) {
            stop("The number of samples in data.extra and in A matrix
                 should be the same!")
        }
        if(length(ex.col) == 1) {
            ex.col <- rep(ex.col, nrow(data.extra))
        } else if (length(ex.col) != nrow(data.extra)) {
            stop("ex.col needs", nrow(data.extra), "colors!")
        }
        if(length(ex.pch) == 1) {
            ex.pch <- rep(ex.pch, nrow(data.extra))
        } else if (length(ex.pch) != nrow(data.extra)) {
            stop("ex.pch needs", nrow(data.extra), "values!")
        }
    }
    
    triA <- function(p){
        a <- norm((p[,2] - p[,1]), '2')
        b <- norm((p[,3] - p[,1]), '2')
        c <- norm((p[,3] - p[,2]), '2')
        s <- (a + b + c) / 2
        return(sqrt(s * (s - a) * (s - b) * (s - c)))
    }
    
    
    A <- A / rep(colSums(A), 1, each = nrow(A))
    cM <- combn(K, 3)
    areaV <- vector('numeric', dim(cM)[2])
    for(i in 1:dim(cM)[2]){
        areaV[i] <- triA(A[, cM[,i]])
    }
    
    cirM <- cbind(cM[, which.max(areaV)], NA)
    al <- norm(A[, cirM[2, 1]] - A[, cirM[1, 1]], '2')
    bl <- norm(A[, cirM[3, 1]] - A[, cirM[2, 1]], '2')
    cl <- norm(A[, cirM[1, 1]] - A[, cirM[3, 1]], '2')
    
    cirM[1, 2] <- 0
    cirM[2, 2] <- al / (al + bl + cl) * 2 * pi
    cirM[3, 2] <- (al + bl) / (al + bl + cl) * 2 * pi
    
    pointInd <- matrix(0, K - 3, 2)
    pointInd[, 1] <- setdiff(1:K, cirM[,1])
    
    for (i in 1:dim(pointInd)[1]){
        for (j in 1:dim(cirM)[1]){
            pointInd[i, 2] <- pointInd[i, 2] + norm(A[, cirM[j, 1]] - A[, pointInd[i, 1]], '2')
        } 
    }
    pointInd <- pointInd[sort.int(pointInd[,2], decreasing = TRUE, index.return = TRUE)$ix, , drop = FALSE]
    
    for (i in 1:dim(pointInd)[1]){
        tempDisM <- matrix(NA, dim(cirM)[1], 2)
        for (j in 1:(dim(tempDisM)[1]-1)){
            tempDisM[j, 1] <- norm(A[, cirM[j, 1]] - A[, pointInd[i, 1]], '2') 
            tempDisM[j, 2] <- norm(A[, cirM[j + 1, 1]] - A[, pointInd[i, 1]], '2')
        }
        tempDisM[dim(tempDisM)[1], 1] <- norm(A[, cirM[1, 1]] - A[, pointInd[i, 1]], '2')
        tempDisM[dim(tempDisM)[1], 2] <- norm(A[, cirM[dim(tempDisM)[1], 1]] - A[, pointInd[i, 1]], '2')
        tInd <- which.min(rowSums(tempDisM))
        if (tInd != dim(tempDisM)[1]){
            tPos <- tempDisM[tInd, 1] * (cirM[tInd + 1, 2] - cirM[tInd, 2]) / 
                (tempDisM[tInd, 1] + tempDisM[tInd, 2]) + cirM[tInd, 2]
            cirM <- rbind(cirM[1:tInd, ], c(pointInd[i, 1], tPos), cirM[(tInd + 1):dim(cirM)[1], ])
        }
        else{
            tPos <- tempDisM[tInd][1] * (2 * pi - cirM[tInd, 2]) / 
                (tempDisM[tInd, 1] + tempDisM[tInd, 2]) + cirM[tInd, 2]
            cirM <- rbind(cirM, c(pointInd[i, 1], tPos))
        }
    }
    tempV <- vector('numeric', dim(cirM)[1])
    for (i in 1:(length(tempV) - 1)){
        tempV[i] <- norm(A[, cirM[i + 1, 1]] - A[, cirM[i, 1]], '2') 
    }
    tempV[length(tempV)] <- norm(A[, cirM[length(tempV), 1]] - A[, cirM[i, 1]], '2')
    suml <- sum(tempV)
    tempV2 <- tempV
    tempV2[1] <- 0
    for (i in 1:(length(tempV) - 1)){
        tempV2[i + 1] <- tempV2[i] + tempV[i] / suml * 2 * pi
    }
    cirM <- cbind(cirM, tempV2)
    
    
    
    PS <- rbind(cos(cirM[, 3]), sin(cirM[, 3]))
    tmp <- PS %*% corpcor::pseudoinverse(A[,cirM[, 1]])
    Xproj <- t(data / rowSums(data))
    Xp <- tmp %*% Xproj
    

    
    plot(Xp[1,], Xp[2,], col = col, cex = cex, pch = pch,
        xlab = NA, ylab = NA, asp = 1, ...)

    for(i in seq_along(MGlist)){
        points(Xp[1,MGlist[[i]]], Xp[2,MGlist[[i]]],
            cex=mg.cex, col=mg.col[i], pch = mg.pch[i] )
    }
    
    # points(PS[1,], PS[2,], pch = 17, xlab = NA, ylab = NA, asp = 1)

    if(!is.null(data.extra)){
        data.extra <- t(data.extra / rowSums(data.extra))
        Xpe <- tmp %*% data.extra
        for(i in seq_len(ncol(data.extra))){
            points(Xpe[1,i], Xpe[2,i],
                   cex=ex.cex, col=ex.col[i], pch = ex.pch[i])
        }
    }
}
