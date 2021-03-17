sffsnnlsFast <- function(Xall, Aall, Kmax = NULL, deltaK = NULL, iInd = NULL){
    Kmax <- max(ncol(Aall), Kmax)
    deltaK <- min(deltaK, ncol(Aall) - Kmax)

    costfunc <- function(combset) {
        A <- Aall[,combset,drop=FALSE]
        fit.err <- apply(Xall, 2, function(x) (nnls::nnls(A,x))$deviance)
        sum(fit.err)
    }
    
    xres <- function(combset) {
        A <- Aall[,combset,drop=FALSE]
        Xres <- apply(Xall, 2, function(x) (nnls::nnls(A,x))$residuals)
        return(Xres)
    }

    adbest <- function(combset, X) {
        if (length(combset) == 1)
            return(combset[1])
        A <- Aall[,combset,drop=FALSE]
        coeff <- apply(X, 2, function(x) (nnls::nnls(A,x))$x)
        coefSum <- sapply(seq_len(dim(A)[2]), function(x) 
            sum((A[, x, drop = FALSE] %*% coeff[x, , drop = FALSE])^2))
        
        
        # show(coefSum)
        return(combset[which.max(coefSum)])
    }
    
    rmworst <- function(combset) {
        A <- Aall[,combset,drop=FALSE]
        coefX <- apply(Xall, 2, function(x) (nnls::nnls(A,x)$x))
        resiA <- sapply(seq_len(dim(A)[2]), function(x) 
            nnls::nnls(A[, -x, drop = FALSE], A[, x, drop = FALSE])$residuals)
        resiX <- Xall - A %*% coefX
        
        deviAll <- sapply(seq_len(dim(A)[2]), function(x)
            sum((resiX + resiA[, x, drop = FALSE] %*%  coefX[x, , drop = FALSE])^2))
        return(combset[which.min(deviAll)])
    }

    allset <- seq_len(ncol(Aall))
    sffsset <- rep(list(NA), Kmax + deltaK)
    sffscost <- rep(Inf, Kmax + deltaK)
    convset <-  c(iInd)
    
    if (length(convset) > 0){
        sffsset[[length(convset)]] <- convset
        sffscost[length(convset)] <- costfunc(convset)
    }

    while(length(convset) < Kmax + deltaK){
        message(convset)
        outset <- setdiff(allset, convset)
        if (length(outset) == 0) break

        #Step 1 (Inclusion)
        
        Xres <- xres(convset)
        addelem <- adbest(outset, Xres)
        convset <- c(convset, addelem)
        newcost <- costfunc(convset)

        if(length(convset) > 2){
            #Step 2 (Conditional exclusion)
            remelem <- rmworst(convset)
            if(addelem != remelem){
                fit.err <- costfunc(setdiff(convset, remelem))
                while(fit.err < sffscost[length(convset)-1]){
                    convset <- setdiff(convset, remelem)
                    newcost <- fit.err
                    if(length(convset) <= 2) break
                    remelem <- rmworst(convset)
                    fit.err <- costfunc(setdiff(convset, remelem))
                }
            }
        }

        sffsset[[length(convset)]] <- convset
        sffscost[length(convset)] <- newcost
    }

    return(list(sffsset=sffsset, sffscost=sffscost))
}
