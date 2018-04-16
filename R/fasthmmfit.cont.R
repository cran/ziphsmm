
#######################################################
#' Fast gradient descent algorithm to learn the parameters
#' in a specialized continuous-time zero-inflated hidden Markov model, where 
#' zero-inflation only happens in State 1. And if there were covariates, they 
#' could only be the same ones for the state-dependent log Poisson means and 
#' the logit structural zero proportion.
#'
#' @param y observed time series values
#' @param x matrix of covariates for the log poisson means and logit zero proportion.
#' Default to NULL.
#' @param M number of latent states
#' @param prior_init a vector of initial values for prior probability for each state
#' @param tpm_init a matrix of initial values for transition rate matrix
#' @param emit_init a vector of initial values for the means for each poisson distribution
#' @param zero_init a scalar initial value for the structural zero proportion 
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
#' tpmparm <- c(-1,-2)
#' zeroindex <- c(1,0)
#' zeroparm <- c(0,-1,1)
#' emitparm <- c(2,0.5,-0.5,3,0.3,-0.2)
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
#'   zero_init=0.4,emit_init=c(7,21), timeindex=timeindex,
#'   hessian=FALSE, method="BFGS", control=list(trace=1))
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to fit a homogeneous hidden semi-markov model
fasthmmfit.cont <- function(y, x=NULL, M, prior_init, tpm_init, emit_init, 
                       zero_init,  yceil=NULL, timeindex,
                       method="Nelder-Mead", hessian=FALSE, ...){
  ntimes <- length(y)
  vdiff <- diff(timeindex)
  udiff <- sort(unique(vdiff))

  if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
  
  if(is.null(x)){
    allparm <- rep(NA, M*M+M)
    allparm[1:(M-1)] <- glogit(prior_init)
    lastindex <- M - 1
    for(i in 1:M){
      for(j in 1:M){
        if(i!=j){
          allparm[lastindex+1] <- glogit(tpm_init[i,j])
          #allparm[lastindex+1] <- log(tpm_init[i,j])
          lastindex <- lastindex + 1
        }
      }
    }
    allparm[lastindex+1] <- log(zero_init) - log(1-zero_init)
    lastindex <- lastindex + 1
    allparm[(lastindex+1):(lastindex+M)] <- log(emit_init)
    
    result <- optim(allparm,zipnegloglik_nocov_cont,grad_zipnegloglik_nocov_cont,
                      M=M,y=y,ntimes=ntimes,timeindex=timeindex, udiff=udiff,
                      method=method,hessian=hessian,...)
      
      fitparm <- result$par
      negloglik <- result$value
      natparm <- retrieve_nocov_cont(fitparm,M)
      return(list(prior=natparm$delta,tpm=natparm$gamma,
                  zeroprop=natparm$theta,emit=natparm$lambda,
                  working_parm=fitparm,negloglik=negloglik))
   
  }else{
    #with covariates
    x <- cbind(1,x)
    ncolx <- ncol(x)
    allparm <- rep(NA, M-1+M*(M-1)+ncolx*(1+M))
    allparm[1:(M-1)] <- glogit(prior_init)
    lastindex <- M - 1
    for(i in 1:M){
      for(j in 1:M){
        if(i!=j){
          allparm[lastindex+1] <- glogit(tpm_init[i,j])
          lastindex <- lastindex + 1
        }
      }
    }
    
    allparm[(lastindex+1):(lastindex+ncolx)] <- 
      c(log(zero_init)-log(1-zero_init),rep(0,ncolx-1))
    lastindex <- lastindex + ncolx
    
    for(i in 1:M){
      allparm[(lastindex+1):(lastindex+ncolx)] <- 
        c(log(emit_init[i]),rep(0,ncolx-1))
      lastindex <- lastindex + ncolx
    }
    
      result <- optim(allparm,zipnegloglik_cov_cont,grad_zipnegloglik_cov_cont,
                      M=M,y=y,covariates=x,ntimes=ntimes,
                      timeindex=timeindex,udiff=udiff,
                      method=method,hessian=hessian,...)
      
      fitparm <- result$par
      negloglik <- result$value
      natparm <- retrieve_cov_cont(fitparm,M,ncolx)
      return(list(prior=natparm$delta,tpm=natparm$gamma,
                  zeroparm=natparm$thetaparm,emitparm=natparm$lambdaparm,
                  working_parm=fitparm,negloglik=negloglik))
  }
}
