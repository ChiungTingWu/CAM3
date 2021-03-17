sffsnnls <- function(Xall, Aall, Kmax = NULL, deltaK = NULL, iInd = NULL){
    Kmax <- max(ncol(Aall), Kmax)
    deltaK <- min(deltaK, ncol(Aall) - Kmax)

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
        fit.err<- unlist(lapply(outset, function(x) costfunc(c(convset,x))))
        newelem <- outset[which.min(fit.err)]
        convset <- c(convset, newelem)
        newcost <- min(fit.err)

        if(length(convset) > 2){
            #Step 2 (Conditional exclusion)
            rmresult <- rmworst(convset)
            if(newelem != rmresult$worstelem){
                while(rmresult$cost < sffscost[length(convset)-1]){
                    convset <- setdiff(convset, rmresult$worstelem)
                    newcost <- rmresult$cost
                    if(length(convset) <= 2) break
                    rmresult <- rmworst(convset)
                }
            }
        }

        sffsset[[length(convset)]] <- convset
        sffscost[length(convset)] <- newcost
    }

    return(list(sffsset=sffsset, sffscost=sffscost))
}
