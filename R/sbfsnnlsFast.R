sbfsnnlsFast <- function(Xall, Aall){
    minK <- 1
    
    costfunc <- function(combset) {
        A <- Aall[,combset,drop=FALSE]
        return(sum(apply(Xall, 2, function(x) (nnls::nnls(A,x))$deviance)))
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
    sbfsset <- rep(list(NA), ncol(Aall))
    sbfscost <- rep(Inf, ncol(Aall))
    currset <-  seq_len(ncol(Aall))
    
    sbfsset[[length(currset)]] <- currset
    sbfscost[length(currset)] <- costfunc(currset)
    
    while(length(currset) > minK){
        message(currset)
        #Step 1 (Exclusion)
        remelem <- rmworst(currset)
        currset <- setdiff(currset, remelem)
        outset <- setdiff(allset, currset)
        
        Xres <- xres(currset)
        newcost <- sum(Xres^2)
        # show('werwer')
        
        #Step 2 (Conditional inclusion)
        
        addelem <- adbest(outset, Xres)
        if(remelem != addelem){
            # show(c('diff', addelem))
            fit.err <- costfunc(c(currset, addelem))
            while(!is.null(fit.err) && fit.err < sbfscost[length(currset)+1]){
                # show('add one back')
                currset <- c(currset, addelem)
                outset <- setdiff(allset, currset)
                newcost <- fit.err
                
                Xres <- xres(currset)
                addelem <- adbest(outset, Xres)
                
                fit.err<- costfunc(c(currset, addelem))
                
            }        
            
        }
        
        sbfsset[[length(currset)]] <- currset
        sbfscost[length(currset)] <- newcost
    }
    
    return(list(sbfsset=sbfsset, sbfscost=sbfscost))
}
