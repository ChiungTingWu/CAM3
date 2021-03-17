#' MG cluster detection for CAM3
#'
#' This function finds corner clusters as MG clusters
#' (clusters containing marker genes).
#' @param PrepResult An object of class "\code{\link[debCAM]{CAMPrepObj}}" 
#'     obtained from \code{\link[debCAM]{CAMPrep}} function.
#' @param fast.mode Use fast mode of greedy search or not. The normal mode may
#'     give more accurate results, but computation time is much longer. The
#'     default is TRUE.
#' @details This function is used internally by \code{\link{CAM3Run}} function
#' to detect clusters containing marker genes,
#' or used when you want to perform CAM step by step.
#'
#' This function provides two solutions. The forward and backward greedy search. 
#' Then the best one of certain source number is selected based on
#' reconstruction errors of all data points in original space.
#' @return A list of class "\code{\link[debCAM]{CAMMGObj}}" containing the
#' following components:
#' \item{idx}{No use in CAM3, will be deprecated in future version.}
#' \item{corner}{The indexes of clusters as detected corners. The first row is
#'     detected by the forward greedy search, and the second row is detected by
#'     the backward greedy search.}
#' \item{error}{Two values. The first one is the reconstruction error of forward
#'     greedy search, and the second one is for the backward greedy search.}
#' @export
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #preprocess data
#' rPrep3 <- CAM3Prep(data, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
#'
#' #Marker gene cluster detection with greedy search
#' rMGC3 <- CAM3MGCluster(rPrep)
CAM3MGCluster <- function(PrepResult, fast.mode=TRUE) {
    X <- PrepResult@centers
  
    Xall <- PrepResult@Xprep[,!(PrepResult@cluster$cluster %in%
                                    PrepResult@c.outlier)]
   
    MG2C <- debCAM::CAMMGCluster(2, PrepResult)
    iInd <- c(which(MG2C@corner[2, 1] == colnames(PrepResult@centers)),
            which(MG2C@corner[2, 2] == colnames(PrepResult@centers)))
    
    if (fast.mode){
        reconsErrB <- sbfsnnlsFast(Xall, X)
        reconsErrF <- sffsnnlsFast(Xall, X, iInd = iInd)
    } else {
        reconsErrB <- sbfsnnls(Xall, X)
        reconsErrF <- sffsnnls(Xall, X, iInd = iInd)
    }
    reconsErrF$sffsset[[1]] <- reconsErrB$sbfsset[[1]]
    
    MGk.list <- vector("list", ncol(X))
    for (k in seq_along(MGk.list)) {
        MGk.list[[k]] <- new("CAMMGObj", idx=c(1, 1), 
            corner = rbind(as.integer(colnames(X[,reconsErrF$sffsset[[k]]])),
                         as.integer(colnames(X[,reconsErrB$sbfsset[[k]]]))),
            error = cbind(reconsErrF$sffscost[k], reconsErrB$sbfscost[k]))
    }

    return(MGk.list)
}