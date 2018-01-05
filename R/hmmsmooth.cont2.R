
########################################################
#' Compute the posterior state probabilities for continuous-time
#' hidden Markov models where zero-inflation only happens in state 1 and
#' covariates can only be included in the state-dependent parameters
#' @param y the observed series to be decoded
#' @param x matrix of covariates for the log poisson means and logit zero proportion.
#' Default to NULL.
#' @param M number of latent states
#' @param prior prior parameters from the fitted continuous-time hidden Markov model
#' @param tpm transition rate parameters from the fitted continuous-time hidden Markov model
#' @param zeroparm parameters for the structural zero proportions in the fitted
#' continuous-time hidden Markov model
#' @param emitparm parameters for the Poisson means in the fitted continuous-time
#' hidden Markov model
#' @param timeindex a vector containing the time points
#' @return posterior state probabilities
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' set.seed(2910)
#' priorparm <- 0
#' tpmparm <- c(-1,-2)
#' zeroindex <- c(1,0)
#' zeroparm <- c(0,-1,1)
#' emitparm <- c(4,0.2,-0.2,5,0.1,-0.1)
#' workparm <- c(priorparm,tpmparm,zeroparm,emitparm)
#' timeindex <- rep(1,1000)
#' for(i in 2:1000) timeindex[i] <- timeindex[i-1] + sample(1:4,1)
#'
#' designx <- matrix(rnorm(2000),nrow=1000,ncol=2)

#' result <- hmmsim2.cont(workparm,2,1000,zeroindex,emit_x=designx,
#'                       zeroinfl_x=designx,timeindex=timeindex)
#' y <- result$series
#' state <- result$state
#'
#' fit2 <-  fasthmmfit.cont(y=y,x=designx,M=2,prior_init=c(0.5,0.5),
#'   tpm_init=matrix(c(-0.2,0.2,0.1,-0.1),2,2,byrow=TRUE),
#'   zero_init=0.4,emit_init=c(50,150), timeindex=timeindex,
#'   hessian=FALSE, method="BFGS", control=list(trace=1))

#' post <- hmmsmooth.cont2(y,designx,2,fit2$prior,fit2$tpm,fit2$zeroparm,
#'        fit2$emitparm,timeindex)
#' 
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hmmsmooth.cont2 <- function(y, x, M, prior, tpm, 
                           zeroparm, emitparm,timeindex){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")

  ntimes <- length(y)
  vdiff <- diff(timeindex)
  udiff <- sort(unique(vdiff))
  expms <- getallexpm(tpm, udiff)
  x <- cbind(1,x)
  fb <- smooth_cov_cont(prior,tpm,zeroparm,emitparm,
                          y,x,ntimes,timeindex,udiff,expms)
  
  return(list(prob_each=fb$Gamma,prob_sum=fb$colsumgamma))
}