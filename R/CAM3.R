#' CAM3: A package for fully unsupervised deconvolution of complex tissues; a
#' new and fast version of debCAM
#'
#' The core function in this package is \code{\link{CAM3Run}} which achieves 
#' fully unsupervised deconvolution on mixture expression profiles.
#' Each step in \code{\link{CAM3Run}} can also be performed separately by
#' \code{\link{CAM3Prep}}, \code{\link{CAM3MGCluster}}
#' and \code{\link{CAM3ASest}} in a more flexible workflow.
#' \code{\link{cotMG}} can help extract a complete marker list from CAM
#' results. \code{\link[debCAM]{MDL}} can help decide the underlying subpopulation
#' number.
#'
#' @docType package
#' @name CAM3-package
#' @aliases CAM3
#'
#' @import methods
#' @import debCAM
#' @import stats
#' @import graphics
#' @importFrom NMF .fcnnls
#' @importFrom corpcor pseudoinverse
NULL
