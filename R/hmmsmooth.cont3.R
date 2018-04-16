
########################################################
#' Compute the posterior state probabilities for continuous-time
#' hidden Markov models with
#' covariates in the state-dependent parameters and transition rates
#' @param y the observed series to be decoded
#' @param x matrix of covariates in the state-dependent parameters and transition rates.
#' @param M number of latent states
#' @param prior prior parameters from the fitted continuous-time hidden Markov model
#' @param tpmparm parameters from the fitted continuous-time hidden Markov model
#' @param zeroparm parameters for the structural zero proportions in the fitted
#' continuous-time hidden Markov model
#' @param emitparm parameters for the Poisson means in the fitted continuous-time
#' hidden Markov model
#' @param timeindex a vector containing the time points
#' @return posterior state probabilities
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' \dontrun{
#' set.seed(2910)
#' priorparm <- 0
#' tpmparm <- c(-2,0.1,-0.1,-2,-0.2,0.2)
#' zeroindex <- c(1,0)
#' zeroparm <- c(0,-1,1)
#' emitparm <- c(2,0.5,-0.5,3,0.3,-0.2)
#' workparm <- c(priorparm,tpmparm,zeroparm,emitparm)
#' timeindex <- rep(1,1000)
#' for(i in 2:1000) timeindex[i] <- timeindex[i-1] + sample(1:4,1)
#'
#' designx <- matrix(rnorm(2000),nrow=1000,ncol=2)
#' result <- hmmsim3.cont(workparm,2,1000,zeroindex,x=designx,timeindex=timeindex)
#' y <- result$series
#' state <- result$state
#'
#' fit2 <-  fasthmmfit.cont3(y=y,x=designx,M=2,
#'   initparm=workparm, timeindex=timeindex,
#'   hessian=FALSE, method="CG", control=list(trace=1))

#' post <- hmmsmooth.cont3(y,designx,2,fit2$prior,fit2$tpm,fit2$zeroparm,
#'        fit2$emitparm,timeindex)
#' }
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hmmsmooth.cont3 <- function(y, x, M, prior, tpmparm, 
                            zeroparm, emitparm,timeindex){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  
  ntimes <- length(y)

  x <- cbind(1,x)
  
  cube <- getallexpm3(M,tpmparm,x,timeindex)
  
  nodeprob <- getnodeprob_cov_cont(y, x, zeroparm, emitparm, M)
  
  fb <- fb_cont3(prior, cube, nodeprob,
                 ntimes, ntimes,timeindex)
  return(list(prob_each=fb$Gamma,prob_sum=fb$colsumgamma))
}