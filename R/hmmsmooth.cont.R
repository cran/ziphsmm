
########################################################
#' Compute the posterior state probabilities for continuous-time
#' hidden Markov models without covariates where zero-inflation only
#' happens in state 1
#' @param y the observed series to be decoded
#' @param M number of latent states
#' @param prior_init a vector of prior probability values
#' @param tpm_init transition rate matrix
#' @param emit_init a vector containing means for each poisson distribution
#' @param zero_init a scalar containing structural zero proportion in state 1
#' @param timeindex a vector containing the time points
#' @return posterior state probabilities
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(10,40,70)
#' zero_init <- c(0.5,0,0)
#' omega <- matrix(c(-0.3,0.2,0.1,0.1,-0.2,0.1,0.2,0.2,-0.4),3,3,byrow=TRUE)
#' timeindex <- rep(1,1000)
#' for(i in 2:1000) timeindex[i] <- timeindex[i-1] + sample(1:3,1)
#' result <- hmmsim.cont(n=1000,M=3,prior=prior_init, tpm_parm=omega,
#'           emit_parm=emit_init,zeroprop=zero_init,timeindex=timeindex)
#' y <- result$series
#' fit2 <-  fasthmmfit.cont(y,x=NULL,M=3,prior_init,omega,
#'                         emit_init,0.5,timeindex=timeindex,hessian=FALSE,
#'                         method="BFGS", control=list(maxit=500,trace=1))
#' post <- hmmsmooth.cont(y,3,fit2$prior,fit2$tpm,fit2$emit,
#'            fit2$zeroprop,timeindex=timeindex)
#' 
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hmmsmooth.cont <- function(y, M, prior_init, tpm_init, 
                            emit_init, zero_init, timeindex){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(prior_init)!=M | length(emit_init)!=M | 
     nrow(tpm_init)!= M | ncol(tpm_init)!=M) stop("The dimension of the initial value does not equal M!") 
  
  ntimes <- length(y)
  vdiff <- diff(timeindex)
  udiff <- sort(unique(vdiff))
  expms <- getallexpm(tpm_init, udiff)
  
  fb <- smooth_nocov_cont(prior_init,tpm_init,zero_init,emit_init,
                          y,ntimes,timeindex,udiff,expms)
  
  return(list(prob_each=fb$Gamma,prob_sum=fb$colsumgamma))
}