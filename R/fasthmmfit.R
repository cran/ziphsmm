
#######################################################
#' Fast gradient descent / stochastic gradient descent algorithm to learn the parameters
#' in a specialized zero-inflated hidden Markov model, where zero-inflation only happens
#' in State 1. And if there were covariates, they could only be the same ones for the
#' state-dependent log Poisson means and the logit structural zero proportion.
#'
#' @param y observed time series values
#' @param x matrix of covariates for the log poisson means and logit zero proportion.
#' Default to NULL.
#' @param ntimes a vector specifying the lengths of individual, 
#' i.e. independent, time series. If not specified, the responses are assumed to 
#' form a single time series, i.e. ntimes=length(y).
#' @param M number of latent states
#' @param prior_init a vector of initial values for prior probability for each state
#' @param tpm_init a matrix of initial values for transition probability matrix
#' @param emit_init a vector of initial values for the means for each poisson distribution
#' @param zero_init a scalar initial value for the structural zero proportion 
#' @param yceil a scalar defining the ceiling of y, above which the values will be
#' truncated. Default to NULL.
#' @param stochastic Logical. Should the stochastic gradient descent methods be used.
#' @param nmin a scalar for the minimum number of observations before the first
#' iteration of stochastic gradient descent. Default to 1000.
#' @param nupdate a scalar specifying the total number of updates for stochastic
#' gradient descent. Default to 100.
#' @param power a scalar representing the power of the learning rate, which should lie
#' between (0.5,1]. Default to 0.7
#' @param rate a vector of learning rate in stochastic gradient descent for the logit
#' parameters and log parameters. Default to c(1,0.05).
#' @param method method to be used for direct numeric optimization. See details in
#' the help page for optim() function. Default to Nelder-Mead.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' Note that the hessian is for the working parameters, which are the generalized logit of prior 
#' probabilities (except for state 1), the generalized logit of the transition probability 
#' matrix(except 1st column), the logit of non-zero zero proportions, and the log of 
#' each state-dependent poisson means
#' @param ... Further arguments passed on to the optimization methods
#' @return the maximum likelihood estimates of the zero-inflated hidden Markov model
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' #1. no covariates
#' set.seed(135)
#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(10,40,70)
#' zero_init <- c(0.5,0,0)
#' omega <- matrix(c(0.5,0.3,0.2,0.4,0.3,0.3,0.2,0.4,0.4),3,3,byrow=TRUE)
#' result <- hmmsim(n=10000,M=3,prior=prior_init, tpm_parm=omega,
#'                 emit_parm=emit_init,zeroprop=zero_init)
#' y <- result$series
#' 
#' time <- proc.time()
#' fit1 <-  fasthmmfit(y,x=NULL,ntimes=NULL,M=3,prior_init,omega,
#'               emit_init,0.5, hessian=FALSE,
#'               method="BFGS", control=list(trace=1))
#' proc.time() - time
#'
#' 
#' #2. with covariates
#' priorparm <- 0
#' tpmparm <- c(-2,2)
#' zeroindex <- c(1,0)
#' zeroparm <- c(0,-1,1)
#' emitparm <- c(2,0.5,-0.5,3,0.3,-0.2)
#' workparm <- c(priorparm,tpmparm,zeroparm,emitparm)
#'
#' designx <- matrix(rnorm(20000),nrow=10000,ncol=2)
#' x <- cbind(1,designx) #has to make the additional 1st column of 1 for intercept
#' result <- hmmsim2(workparm,2,10000,zeroindex,emit_x=designx,zeroinfl_x=designx)
#' y <- result$series
#' 
#' time <- proc.time()
#' fit2 <-  fasthmmfit(y=y,x=x,ntimes=NULL,M=2,prior_init=c(0.5,0.5),
#'               tpm_init=matrix(c(0.9,0.1,0.1,0.9),2,2),
#'               zero_init=0.4,emit_init=c(7,21), hessian=FALSE,
#'               method="BFGS", control=list(trace=1))
#' proc.time() - time
#' fit2
#' 
#' #3. stochastic gradient descent without covariates
#' #no covariates
#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(10,40,70)
#' zero_init <- c(0.5,0,0)
#' omega <- matrix(c(0.5,0.3,0.2,0.4,0.3,0.3,0.2,0.4,0.4),3,3,byrow=TRUE)
#' result <- hmmsim(n=50000,M=3,prior=prior_init, tpm_parm=omega,
#'                 emit_parm=emit_init,zeroprop=zero_init)
#' y <- result$series
#'
#' initparm2 <- c(-1,-0.5,  -0.3,-0.3,-0.4,-0.4,0.5,0.5,  0,2,3,4)
#' time <- proc.time()
#' fitst <- fasthmmfit(y=y,x=NULL,ntimes=NULL,M=3,prior_init=c(0.4,0.3,0.3),
#'               tpm_init=matrix(c(0.6,0.3,0.1,0.3,0.4,0.3,0.1,0.3,0.6),3,3,byrow=TRUE),
#'               zero_init=0.3,emit_init=c(8,35,67),stochastic=TRUE,
#'               nmin=1000,nupdate=1000,power=0.6,rate=c(1,0.05))
#' proc.time() - time
#' str(fitst)
#' 
#' #with covariates
#' priorparm <- 0
#' tpmparm <- c(-2,2)
#' zeroindex <- c(1,0)
#' zeroparm <- c(0,-1,1)
#' emitparm <- c(2,0.5,-0.5,3,0.3,-0.2)
#' workparm <- c(priorparm,tpmparm,zeroparm,emitparm)
#'
#' designx <- matrix(rnorm(100000),nrow=50000,ncol=2)
#' x <- cbind(1,designx) #has to make the additional 1st column of 1 for intercept
#' result <- hmmsim2(workparm,2,50000,zeroindex,emit_x=designx,zeroinfl_x=designx)
#' y <- result$series
#'
#' initparm <- c(0, -1.8,1.8, 0,-0.8,0.8, 1.8,0.6,-0.6,3.1,0.4,-0.3)
#'
#' time <- proc.time()
#' fitst <- fasthmmfit(y=y,x=x,ntimes=NULL,M=2,prior_init=c(0.4,0.6),
#'               tpm_init=matrix(c(0.8,0.2,0.2,0.8),2,2,byrow=TRUE),
#'               zero_init=0.3,emit_init=c(10,25),stochastic=TRUE,
#'               nmin=1000,nupdate=1000,power=0.6,rate=c(1,0.05))
#' proc.time() - time
#' str(fitst)
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to fit a homogeneous hidden semi-markov model
fasthmmfit <- function(y, x=NULL, ntimes=NULL, M, prior_init, tpm_init, emit_init, 
                  zero_init,  yceil=NULL,
                  stochastic=FALSE, nmin=1000, nupdate=100, power=0.7, 
                  rate=c(1,0.05),
                  method="Nelder-Mead", hessian=FALSE, ...){
    if(is.null(ntimes)) ntimes <- length(y)
    if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
      
    if(is.null(x)){
      allparm <- rep(NA, M*M+M)
      allparm[1:(M-1)] <- glogit(prior_init)
      lastindex <- M - 1
      for(i in 1:M){
        allparm[(lastindex + 1):(lastindex + M-1)] <- glogit(tpm_init[i,])
        lastindex <- lastindex + M - 1
      }
      allparm[lastindex+1] <- log(zero_init) - log(1-zero_init)
      lastindex <- lastindex + 1
      allparm[(lastindex+1):(lastindex+M)] <- log(emit_init)
      
      if(stochastic==FALSE){
        result <- optim(allparm,zipnegloglik_nocov,grad_zipnegloglik_nocov,
            M=M,y=y,ntimes=ntimes,
            method=method,hessian=hessian,...)
      
      fitparm <- result$par
      negloglik <- result$value
      natparm <- retrieve_nocov(fitparm,M)
      return(list(prior=natparm$delta,tpm=natparm$gamma,
                  zeroprop=natparm$theta,emit=natparm$lambda,
                  working_parm=fitparm,negloglik=negloglik))
         }else{
           
            ratevec <- c(rep(rate[1],M*M),rep(rate[2],M))
            fitparm <- sgd_zip_nocov(allparm,y,M,ntimes,nmin,nupdate,power,
                                     ratevec)
            natparm <- retrieve_nocov(fitparm[nupdate,],M)
            return(list(prior=natparm$delta,tpm=natparm$gamma,
                        zeroprop=natparm$theta,emit=natparm$lambda,
                        working_parm=fitparm))
         }
    }else{
      x <- cbind(1,x)
      ncolx <- ncol(x)
      allparm <- rep(NA, M-1+M*(M-1)+ncolx*(1+M))
      allparm[1:(M-1)] <- glogit(prior_init)
      lastindex <- M - 1
      for(i in 1:M){
        allparm[(lastindex + 1):(lastindex + M-1)] <- glogit(tpm_init[i,])
        lastindex <- lastindex + M - 1
      }
      
      allparm[(lastindex+1):(lastindex+ncolx)] <- 
        c(log(zero_init)-log(1-zero_init),rep(0,ncolx-1))
      lastindex <- lastindex + ncolx
      
      for(i in 1:M){
        allparm[(lastindex+1):(lastindex+ncolx)] <- 
          c(log(emit_init[i]),rep(0,ncolx-1))
        lastindex <- lastindex + ncolx
      }
      
      if(stochastic==FALSE){
        result <- optim(allparm,zipnegloglik_cov,grad_zipnegloglik_cov,
                      M=M,y=y,covariates=x,ntimes=ntimes,
                      method=method,hessian=hessian,...)
      
        fitparm <- result$par
        negloglik <- result$value
        natparm <- retrieve_cov(fitparm,M,ncolx)
        return(list(prior=natparm$delta,tpm=natparm$gamma,
                  zeroparm=natparm$thetaparm,emitparm=natparm$lambdaparm,
                  working_parm=fitparm,negloglik=negloglik))
        }else{
          ratevec <- c(rep(rate[1],M*M-1+ncolx),
                       rep(rate[2],M*ncolx))
          fitparm <- sgd_zip_cov(allparm,y,x,M,ntimes,nmin,nupdate,power,
                                   ratevec)
          natparm <- retrieve_cov(fitparm[nupdate,],M,ncolx)
          return(list(prior=natparm$delta,tpm=natparm$gamma,
                      zeroprop=natparm$thetaparm,emitparm=natparm$lambdaparm,
                      working_parm=fitparm))
        }
    }
}
