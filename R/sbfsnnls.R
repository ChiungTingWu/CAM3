sbfsnnls <- function(Xall, Aall){
    minK <- 1
    
    costfunc <- function(combset) {
        A <- Aall[,combset,drop=FALSE]
        fit.err <- apply(Xall, 2, function(x) (nnls::nnls(A,x))$deviance)
        sum(fit.err)
    }
    
    rmworst <- function(combset) {
        fit.err <- unlist(
            lapply(seq_along(combset), function(x) costfunc(combset[-x]))
        )
        list(worstelem = combset[which.min(fit.err)], cost = min(fit.err))
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
        rmresult <- rmworst(currset)
        remelem <- rmresult$worstelem
        currset <- setdiff(currset, remelem)
        outset <- setdiff(allset, currset)
        newcost <- rmresult$cost
        #fit.err<- unlist(lapply(outset, function(x) costfunc(c(convset,x))))
        #newelem <- outset[which.min(fit.err)]
        #convset <- c(convset, newelem)
        #newcost <- min(fit.err)
        
       
        #Step 2 (Conditional inclusion)
        fit.err<- unlist(lapply(outset, function(x) costfunc(c(currset,x))))
        addelem <- outset[which.min(fit.err)]
        if(remelem != addelem){
            while(!is.null(fit.err) && min(fit.err) < sbfscost[length(currset)+1]){
                currset <- c(currset, addelem)
                outset <- setdiff(allset, currset)
                newcost <- min(fit.err)
                fit.err<- unlist(lapply(outset, function(x) costfunc(c(currset,x))))
                addelem <- outset[which.min(fit.err)]
            }        
            
        }
        
        #rmresult <- rmworst(convset)
        #if(newelem != rmresult$worstelem){
        #    while(rmresult$cost < sffscost[length(convset)-1]){
        #        convset <- setdiff(convset, rmresult$worstelem)
        #        newcost <- rmresult$cost
        #        if(length(convset) <= 2) break
        #        rmresult <- rmworst(convset)
        #    }
        #}

        
        sbfsset[[length(currset)]] <- currset
        sbfscost[length(currset)] <- newcost
    }
    
    return(list(sbfsset=sbfsset, sbfscost=sbfscost))
}