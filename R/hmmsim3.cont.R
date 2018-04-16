
#######################################################
#' Simulate a continuous-time hidden Markov series and its underlying 
#' states with covariates in state-dependent parameters and transition rates.
#' @param workparm working parameters
#' @param M number of latent states
#' @param n length of the simulated series
#' @param zeroindex a vector specifying whether a certain state is zero-inflated
#' @param x matrix of covariates for the log poisson means, structural zero proportions
#' and transition rates. 
#' @param timeindex a vector containing the time points
#' @return simulated series and corresponding states
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' priorparm <- 0
#' tpmparm <- c(-2,0.1,-0.1,-2,0.2,-0.2)
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
#'state <- result$state
#'
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export

hmmsim3.cont <- function(workparm, M, n, zeroindex, x, timeindex){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  
  #check if covariates
  tempmat <- matrix(0,nrow=n,ncol=1) #for NULL covariates
  ncovprior <- 0
  covprior <- tempmat
  
  ncovtpm <- ncol(x)
  ncovemit <- ncol(x)
  ncovzeroinfl <- ncol(x)
  covtpm <- x
  covemit <- x
  covzeroinfl <- x

  
    result <- hmm_cov_gen_cont3(workparm, M, n, ncovprior, covprior, ncovtpm, covtpm,
                               ncovzeroinfl, covzeroinfl, ncovemit, covemit, zeroindex, timeindex)
    series <- result[,1]
    state <- result[,2]
    
    return(list(series=series,state=state))
    
  
}
