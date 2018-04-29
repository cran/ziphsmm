#' \tabular{ll}{
#' Package: \tab ziphsmm \cr
#' Type: \tab Package\cr
#' Version: \tab 2.0.5\cr
#' Date: \tab 2018-04-28\cr
#' License: \tab GPL-2\cr
#' LazyLoad: \tab yes\cr
#' LazyData: \tab yes\cr
#' }
#'
#' @author Zekun Xu \email{zekunxu@gmail.com}
#' @author Ye Liu
#' Maintainer: Zekun Xu \email{zekunxu@gmail.com}
#' @name package-ziphsmm
#' @aliases ziphsmm-package
#' @docType package
#' @title zero-inflated poisson hidden (semi-)Markov models
#' @keywords zero-inflated poisson, hidden Markov models
NULL

##########################################################

#' @useDynLib ziphsmm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics points
#' @importFrom stats optim
#' @importFrom stats runif
#' @importFrom pracma grad

#generalized logit and inverse logit function
#1 / (1+sum(exp(x[-1]))) = p1
#exp(x[k]) / (1+sum(exp(x[-1]))) = pk


glogit <- function(p){
  k <- length(p) - 1
  if(k==0) {x <- log(p) - log(1-p)}else{
  x <- rep(NA, k)
  for(j in 1:k) x[j] <- log(p[j+1]) - log(p[1])}
  return(x)
}


ginvlogit <- function(x){
  k <- length(x) + 1
  p <- rep(NA,k)
  den <- 1+sum(exp(x))
  p[1] <- 1 / den
  for(j in 2:k) p[j] <- exp(x[j-1]) / den
  return(p)
} 






