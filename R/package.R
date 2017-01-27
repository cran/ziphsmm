#' \tabular{ll}{
#' Package: \tab ziphsmm \cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.2\cr
#' Date: \tab 2017-01-26\cr
#' License: \tab GPL-2\cr
#' LazyLoad: \tab yes\cr
#' LazyData: \tab yes\cr
#' }
#'
#' @author Zekun Xu \email{zekunxu@gmail.com}
#' @author Ye Liu \email{yliu87@ncsu.edu}
#' Maintainer: Zekun Xu \email{zekunxu@gmail.com}
#' @name package-ziphsmm
#' @aliases ziphsmm-package
#' @docType package
#' @title zero-inflated poisson hidden (semi-)Markov models
#' @keywords zero-inflated poisson, hidden Markov models
NULL

##########################################################

#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @importFrom graphics points
#' @importFrom stats optim

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






