#' Data preprocessing for CAM3
#'
#' This function perform preprocessing for CAM, including norm-based filtering,
#' dimension deduction, perspective projection, local outlier removal and
#' aggregation of gene expression vectors by clustering.
#' @param data Matrix of mixture expression profiles.
#'     Data frame, SummarizedExperiment or ExpressionSet object will be
#'     internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     Data should be in non-log linear space with non-negative numerical values
#'     (i.e. >= 0). Missing values are not supported.
#'     All-zero rows will be removed internally.
#' @param dim.rdc Reduced data dimension; should be not less than maximum
#'     candidate K.
#' @param thres.low The lower bound of percentage of genes to keep for CAM
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 0.05.
#' @param thres.high The higher bound of percentage of genes to keep for CAM
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 1.
#' @param cluster.method The method to do clustering.
#'     The default "Fixed-Radius" will make all the clusters with the same size.
#'     The alternative "K-Means" will use \code{\link{kmeans}}.
#' @param radius.thres The "cosine" radius of "Fixed-Radius" clustering. The
#'     default is 0.95
#' @param sim.thres The cosine similarity threshold of cluster centers. For
#'     clusters with cosine similarity higher than the threshold, they would be
#'     merged until the number of clusters equals to cluster.num. This parameter
#'     could control the upper bound of similarity amoung sources. The default 
#'     is 0.95.
#' @param cluster.num The lower bound of cluster number, which should be much 
#'     larger than K. The default is 50.
#' @param MG.num.thres The clusters with the gene number smaller than
#'     MG.num.thres will be treated as outliers. The default is 20.
#' @param sample.weight Vector of sample weights. If NULL, all samples have
#'     the same weights. The length should be the same as sample numbers.
#'     All values should be positive.
#' @details This function is used internally by \code{\link{CAM3Run}} function
#' to preprocess data, or used when you want to perform CAM step by step.
#'
#' Low/high-expressed genes are filtered by their L2-norm ranks.
#' Dimension reduction is slightly different from PCA.
#' The first loading vector is forced to be c(1,1,...,1) with unit norm
#' normalization. The remaining are eigenvectors from PCA in the space
#' orthogonal to the first vector.
#' Perspective projection is to project dimension-reduced gene expression
#' vectors to the hyperplane orthogonal to c(1,0,...,0), i.e., the first axis
#' in the new coordinate system.
#' Finally, gene expression vectors are aggregated by clustering
#' to further reduce the impact of noise/outlier and help improve the efficiency
#' of simplex corner detection.
#' @return An object of class "\code{\link[debCAM]{CAMPrepObj}}" containing the
#' following components:
#' \item{Valid}{logical vector to indicate the genes left after filtering.}
#' \item{Xprep}{Preprocessed data matrix.}
#' \item{Xproj}{Preprocessed data matrix after perspective projection.}
#' \item{W}{The matrix whose rows are loading vectors.}
#' \item{SW}{Sample weights.}
#' \item{cluster}{cluster results including two vectors.
#'     The first indicates the cluster to which each gene is allocated.
#'     The second is the number of genes in each cluster.}
#' \item{c.outlier}{The clusters with the gene number smaller than
#'     MG.num.thres.}
#' \item{centers}{The centers of candidate corner clusters (candidate clusters
#'     containing marker genes).}
#' @export
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #preprocess data
#' rPrep3 <- CAM3Prep(data, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
CAM3Prep <- function(data, dim.rdc = 10, thres.low = 0.05, thres.high = 1,
                    cluster.method = c('Fixed-Radius' ,'K-Means'),
                    radius.thres = 0.95, sim.thres = 0.95, cluster.num = 50,
                    MG.num.thres = 20, sample.weight = NULL) {
                    
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
    if (sum(is.na(data)) > 0) {
        stop("Data with missing values are not supported!")
    }
    if (sum(data<0) > 0) {
        stop("Only non-negative data are supported!")
    }
    if (thres.low >= thres.high) {
        stop("thres.low must be smaller than thres.high!")
    }
    if (is.null(rownames(data))) {
        rownames(data) <- seq_len(nrow(data))
    }
    sample.weight <- rep(1, ncol(data))
    data <- data[rowSums(data) > 0,]
    X <- t(data)
    X <- X * sample.weight
    
    ############### Low & high expression filtering ###########################
    sigNorm <- apply(X, 2, function(x) norm(matrix(x),"F") )
    Valid <- sigNorm >= quantile(sigNorm, thres.low) &
        sigNorm <= quantile(sigNorm, thres.high)
    X <- X[,Valid]
    
    Lorg <- dim(X)[1]
    dataSize<- dim(X)[2]
    Xorg <- X
    
    
    ############### PCA dimension deduction ###################################
    L <- min(dim.rdc, Lorg)
    Xmean <- rowMeans(Xorg)
    Xmeanrm <- Xorg - Xmean
    
    #use sample.weight as the first dimension, then PCA
    P1 <- sample.weight
    P1 <- P1/ sqrt(sum(P1^2))
    C1 <- P1 %*% Xmeanrm
    Xmeanrm <- Xmeanrm - matrix(rep(C1, 1, each=Lorg), Lorg) *
        matrix(P1, Lorg, dataSize)
    r <- svd(Xmeanrm)
    Ppca <- t(r$u[,seq_len(L-1)])
    X <- rbind(P1,Ppca) %*% Xorg
    weightMatrix <- rbind(P1, Ppca)
    
    
    ############### perspective projection ####################################
    Xcenter <- c(1, rep(0, L-1))
    denom <- Xcenter %*% X
    denom <- denom[rep(1,L),]
    Xproj <- X / denom
    
    
    ############### local outlier removal #####################################
    XcosSim <- cosSim(Xproj)
 
    Xnei <- XcosSim >= radius.thres
    Xnei.num <- rowSums(Xnei)
    cos.outlier <- Xnei.num < MG.num.thres
    Valid[Valid==TRUE][cos.outlier]<- FALSE
    X <- X[,!cos.outlier]
    Xproj <- Xproj[,!cos.outlier]
    Xnei <- Xnei[!cos.outlier ,!cos.outlier]
    
    
    ############### clustering ################################################
    if (length(cluster.method) > 1) 
        cluster.method <- cluster.method[1]
    
    if (cluster.method == 'Fixed-Radius') {
        clus.list <- list()
        no.clus.ind <- rep(TRUE, dim(Xnei)[1])
        names(no.clus.ind) <- rownames(Xnei) 
        Xnei.num <- rowSums(Xnei[no.clus.ind, no.clus.ind])
        clus.ind <- rep(-1, dim(Xnei)[1])
        count = 1
        while (sum(Xnei.num >= MG.num.thres) > 0){
            newClusN <- which.max(Xnei.num)
            newClus <- which(Xnei[no.clus.ind, names(newClusN)])
            clus.ind[which(no.clus.ind)[newClus]] <- count
            if (length(newClus) == 1)
                newClus <- newClusN
            clus.list[[names(newClusN)]] <- names(newClus)
            no.clus.ind[names(newClus)] <- FALSE
            Xnei.num <- rowSums(Xnei[no.clus.ind, no.clus.ind])
            
            show(count)
            show(c(sum(Xnei.num >= MG.num.thres), sum(no.clus.ind)))
            show(Sys.time())
            count = count + 1
        }
        cluster <- list(cluster = clus.ind, size = sapply(clus.list, length))
    } else {
        if (cluster.method != 'K-Means')
            stop("Only Fixed-Radius and K-Means are supported.")
        
        clusterRes <- kmeans(t(Xproj), cluster.num, iter.max=100)
        #repeat K-means and use the best one
        if (ncol(Xproj) > 15000) {
            ntry <- 10
        } else {
            ntry <- 50
        }
        for (i in seq_len(ntry)){
            tmp <- kmeans(t(Xproj), cluster.num, iter.max=100)
            if (clusterRes$tot.withinss > tmp$tot.withinss){
                clusterRes <- tmp
            }
        }
        cluster <- list(cluster = clusterRes$cluster, size = clusterRes$size)
    }
    
    cluster.valid <- which(cluster$size >= MG.num.thres)
    c.outlier <- which(cluster$size < MG.num.thres)
    message('outlier cluster number: ',sum(cluster$size < MG.num.thres),"\n")

    medcenters <- sapply(cluster.valid, function(x)
        rowMeans(Xproj[, which(cluster$cluster==x)]))
    colnames(medcenters) <- cluster.valid
    
    
    ################ convex hull by linear programming ########################
    message("Linear programming...\n")
    
    obval <- rep(NA, dim(medcenters)[2])
    f.con <- rbind(medcenters, rep(1, dim(medcenters)[2]))
    f.dir <- rep('=', L + 1)
    for(n in 1:dim(medcenters)[2]){
        f.obj <- rep(0, dim(medcenters)[2])
        f.obj[n] <- 1
        f.rhs <- c(medcenters[,n], 1)
        obval[n] <- lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)$objval
    }
    corner <- which(obval > 0.5)
    message('Convex hull cluster number: ',length(corner),"\n")
    
    
    ################ clusters merging #########################################
    
    centers <- medcenters[, corner]
    clus.size <- sapply(clus.list[corner], length)
    clus.ind <- vector("list", ncol(centers))
    names(clus.ind) <- colnames(centers)
    for (x in seq_along(clus.ind))
        clus.ind[[x]] <- as.integer(names(clus.ind[x]))
    
    names(clus.size) <- colnames(centers)
    cen.cos <- cosSim(centers)
    diag(cen.cos) <- NA

    while(max(cen.cos, na.rm=TRUE) > sim.thres && ncol(centers) > cluster.num){
        maxInd <- which(cen.cos == max(cen.cos, na.rm=TRUE), arr.ind=TRUE)[1, ]
        clus1 <- rownames(cen.cos[maxInd[1], maxInd[2], drop = FALSE]) 
        clus2 <- colnames(cen.cos[maxInd[1], maxInd[2], drop = FALSE])
        
        newCenter <- (centers[,clus1, drop=FALSE] * clus.size[clus1] + 
                          centers[,clus2, drop=FALSE] * clus.size[clus2]) / 
            (clus.size[clus1] + clus.size[clus2]) 
        newClusSize <- clus.size[clus1] + clus.size[clus2]
        
        centers <- cbind(centers, newCenter)
        centers <- centers[,-maxInd]
        clus.size <- c(clus.size, newClusSize)
        clus.size <- clus.size[-maxInd]
        new.ind <- list(c(clus.ind[[clus1]], clus.ind[[clus2]]))
        names(new.ind) <- clus1
        clus.ind <- clus.ind[-maxInd]
        clus.ind <- c(clus.ind, new.ind)
        
        cen.cos <- cosSim(centers)
        diag(cen.cos) <- NA
    }
    
    for (x in seq_along(clus.ind)) {
        if (length(clus.ind[[x]]) > 1) {
            for (i in clus.ind[[x]]) {
                show(i)
                cluster$cluster[cluster$cluster == as.integer(i)] <- as.integer(names(clus.ind[x]))
            }
        }
    }
    
    return(new("CAMPrepObj", Valid=Valid, Xprep=X, Xproj=Xproj, W=weightMatrix,
               SW=sample.weight, cluster=cluster, c.outlier=c.outlier,
               centers=centers))
}