
#######################################################
#' Fast gradient descent algorithm to learn the parameters
#' in a specialized continuous-time zero-inflated hidden Markov model, where 
#' zero-inflation only happens in State 1 with covariates in the state-dependent 
#' parameters and transition rates.
#'
#' @param y observed time series values
#' @param x matrix of covariates in the state-dependent 
#' parameters and transition rates.
#' @param M number of latent states
#' @param initparm vector of initial working parameters for prior, transition, 
#' zero proportion, and emission parameters. 
#' @param yceil a scalar defining the ceiling of y, above which the values will be
#' truncated. Default to NULL.
#' @param timeindex a vector containing the time points
#' @param method method to be used for direct numeric optimization. See details in
#' the help page for optim() function. Default to Nelder-Mead.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' Note that the hessian is for the working parameters, which are the generalized logit of prior 
#' probabilities (except for state 1), the generalized logit of the transition probability 
#' matrix(except 1st column), the logit of non-zero zero proportions, and the log of 
#' each state-dependent poisson means
#' @param ... Further arguments passed on to the optimization methods
#' @return the maximum likelihood estimates of the zero-inflated hidden Markov model
#' @references Liu, Yu-Ying, et al. "Efficient learning of continuous-time hidden 
#' markov models for disease progression." Advances in neural information 
#' processing systems. 2015.
#' @examples
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
#' fit2 <-  fasthmmfit.cont3(y=y,x=designx,M=2,initparm=workparm,
#'   timeindex=timeindex,
#'   hessian=FALSE, method="BFGS", control=list(trace=1))
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to fit a homogeneous hidden semi-markov model
fasthmmfit.cont3 <- function(y, x, M, initparm,  yceil=NULL, timeindex,
                            method="Nelder-Mead", hessian=FALSE, ...){
  ntimes <- length(y)
  vdiff <- diff(timeindex)
  udiff <- sort(unique(vdiff))
  
  if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
  
    #with covariates
    x <- cbind(1,x)
    ncolx <- ncol(x)
    
    result <- optim(initparm,zipnegloglik_cov_cont3,
                    M=M,y=y,covariates=x,ntimes=ntimes,
                    timeindex=timeindex,
                    method=method,hessian=hessian,...)
    
    fitparm <- result$par
    negloglik <- result$value
    natparm <- retrieve_cov_cont3(fitparm,M,ncolx)
    return(list(prior=natparm$delta,tpmparm=natparm$gammaparm,
                zeroparm=natparm$thetaparm,emitparm=natparm$lambdaparm,
                working_parm=fitparm,negloglik=negloglik))
  
}
