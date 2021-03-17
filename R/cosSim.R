cosSim <- function (x){
    cosUp <- crossprod(x)
    cosDown <- sqrt(colSums(x^2))
    cosDown <- cosDown %*% t(cosDown)
    return(cosUp / cosDown)
}


