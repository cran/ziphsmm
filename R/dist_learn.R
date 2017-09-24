
#######################################################
#' Distributed learning for a longitudinal continuous-time zero-inflated Poisson
#' hidden Markov model, where zero-inflation only happens in State 1. Assume that
#' the intercepts are different among subjects, a distributed learning algorithm
#' is used to compute the common coefficients of thecovariates for the 
#' state-dependent log Poisson means and the logit structural zero proportion.
#' @param ylist list of observed time series values for each subject
#' @param xlist list of design matrices for each subject. 
#' @param timelist list of time indices
#' @param prior_init a vector of initial values for prior probability for each state
#' @param tpm_init a matrix of initial values for transition rate matrix
#' @param emit_init a vector of initial values for the means for each poisson distribution
#' @param zero_init a scalar initial value for the structural zero proportion 
#' @param yceil a scalar defining the ceiling of y, above which the values will be
#' truncated. Default to NULL. 
#' @param rho tuning parameters in the distributed learning algorithm.
#' @param maxit maximum number iteration. Default to 100.
#' @param tol tolerance in the terms of the relative change in the norm of the
#' common coefficients. Default to 1e-4. 
#' @param print whether to print each iteration. Default to TRUE.
#' @return the maximum likelihood estimates of the zero-inflated hidden Markov model
#' @references Liu, Yu-Ying, et al. "Efficient learning of continuous-time hidden 
#' markov models for disease progression." Advances in neural information 
#' processing systems. 2015.
#' @examples
#' \dontrun{
#' set.seed(930518)
#' ylist <- vector(mode="list",length=20)
#' xlist <- vector(mode="list",length=20)
#' timelist <- vector(mode="list",length=20)

#' for(n in 1:20){
#'  priorparm <- 0 + runif(1,-0.1,0.1)
#'  tpmparm <- c(-1+ runif(1,-0.1,0.1),-2+ runif(1,-0.2,0.2))
#'  zeroindex <- c(1,0)
#'  zeroparm <- c(0+ runif(1,-0.1,0.1),-1,1)
#'  emitparm <- c(2+ runif(1,-0.1,0.1),0.5,-0.5,
#'                3+ runif(1,-0.1,0.1),0.3,-0.2)
#'  workparm <- c(priorparm,tpmparm,zeroparm,emitparm)
#'  timeindex <- rep(1,1440)
#'  for(i in 2:1440) timeindex[i] <- timeindex[i-1] + sample(1:4,1)
#'  timelist[[n]] <- timeindex
#'  
#'  xlist[[n]] <- matrix(rnorm(2880),nrow=1440,ncol=2)
#'  
#'  result <- hmmsim2.cont(workparm,2,1440,zeroindex,emit_x=xlist[[n]],
#'                         zeroinfl_x=xlist[[n]],timeindex=timeindex)
#'  ylist[[n]] <- result$series
#'}
#'
#' prior_init=c(0.5,0.5)
#' tpm_init=matrix(c(-0.2,0.2,0.1,-0.1),2,2,byrow=TRUE)
#' zero_init=0.4
#' emit_init=c(7,21)
#'
#' result <- dist_learn(ylist, xlist, timelist, prior_init, tpm_init, 
#'                     emit_init, zero_init, rho=1, 
#'                     maxit=20, tol=1e-4, print=TRUE)
#'
#' }
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export

dist_learn <- function(ylist, xlist, timelist, prior_init, tpm_init, 
                       emit_init, zero_init, yceil=NULL,
                       rho=1, maxit=100, tol=1e-4, print=TRUE){

  nsubj <- length(ylist)
  M <- ncol(tpm_init)
  
  #retrieve working parameters
  ncolx <- ncol(xlist[[1]]) + 1
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
  
  #fixed indices in working parameters
  indices <- c(1:(M*M-1), seq(M*M, by=ncolx, length=M+1))
    
  #paramters and their gradients
  parm <- matrix(0, nsubj, length(allparm))
  for(i in 1:nsubj) parm[i,] <- allparm 
  
  parmbar <- colMeans(parm)
  olddiff <- parm[,-indices,drop=FALSE]-parmbar[-indices]
  olddiff <- sum(olddiff^2)
  
  grad <- matrix(0, nsubj, length(allparm))
  l <- matrix(0,nsubj,length(allparm)-length(indices))
  #ADMM
  #thetai := argmin(fi(thetai) + t(li)%*%(thetai-thetabar) + 0.5*rho*|thetai-thetabar|^2)
  #thetabar := mean(thetai)
  #li := li + rho(thetai - thetabar)
  newf <- function(initparm,y,covariates,M,ntimes,timeindex,udiff,
                   parmbar,indices,rho,l){
    part1 <- zipnegloglik_cov_cont(initparm,y,covariates,M,ntimes,timeindex,udiff)
    diff <- initparm[-indices] - parmbar[-indices]
    part2 <- t(l[i,])%*%diff
    part3 <- 0.5*rho*t(diff)%*%diff
    return(part1+part2+part3)
  }
  
  newgradf <- function(initparm,y,covariates,M,ntimes,timeindex,udiff,
                       parmbar,indices,rho,l){
    part1 <- grad_zipnegloglik_cov_cont(initparm,y,covariates,M,ntimes,timeindex,udiff)
    part2 <- l[i,]
    part3 <- rho * (initparm[-indices] - parmbar[-indices])
    part1[-indices] <- part1[-indices] + part2 + part3 
    return(part1)
  }
  
  iteration <- 1
  nllk <- 0
  #recursion
  while(iteration<=maxit){
    
    #distributed
    oldlik <- nllk
    tempresult <- lapply(1:nsubj, function(i){
      y <- ylist[[i]]
      if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
      x <- cbind(1,xlist[[i]])
      timeindex <- timelist[[i]]
      ntimes <- length(y)
      
      vdiff <- diff(timeindex)
      udiff <- sort(unique(vdiff))
 
      #l[i,] <- l[i,] + rho*(parm[i,-indices] - parmbar[-indices])
      initparm <- parm[i,]
      
      optim(par=initparm,fn=newf,gr=newgradf,
                      M=M,y=y,covariates=x,ntimes=ntimes,
                      timeindex=timeindex,udiff=udiff, 
                      parmbar=parmbar,indices=indices,rho=rho,l=l,
                      method="CG")
                      #method="L-BFGS-B")
     #newf(initparm,y,x,M,ntimes,timeindex,udiff,parmbar,indices,rho,l)
      #grad(newf,initparm,M=M,y=y,covariates=x,ntimes=ntimes, timeindex=timeindex,udiff=udiff, parmbar=parmbar,indices=indices,rho=rho,l=l)
    
      #newgradf(initparm,y,x,M,ntimes,timeindex,udiff,parmbar,indices,rho,l)
    })
    
    nllk <- sum(sapply(1:nsubj,function(i)tempresult[[i]]$value))
    parm <- t(sapply(1:nsubj,function(i) tempresult[[i]]$par))
    
    parmbar <- colMeans(parm)
    
    for(i in 1:nsubj)
      l[i,] <- l[i,] + rho*(parm[i,-indices] - parmbar[-indices])
    
    newdiff <- parm[,-indices,drop=FALSE]-
           matrix(rep(parmbar[-indices],nsubj),
                  nsubj,length(allparm)-length(indices),byrow=TRUE)
    newdiff <- sum(newdiff^2)
    
    relchange <- newdiff/(1+olddiff)
    #if(iteration==maxit | (iteration>3 & relchange<tol)) {
      #for(nn in 1:nsubj) parm[nn, -indices] <- parmbar[-indices]
      #nllk <- sapply(1:nsubj, function(i){
      #  y <- ylist[[i]]
      #  if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
      #  x <- cbind(1,xlist[[i]])
      #  timeindex <- timelist[[i]]
      #  ntimes <- length(y)
        
      #  vdiff <- diff(timeindex)
      #  udiff <- sort(unique(vdiff))
        
      #  l[i,] <- l[i,] + rho*(parm[i,-indices] - parmbar[-indices])
      #  initparm <- parm[i,]
      #  
      #  zipnegloglik_cov_cont(parm[i,],y,x,M,ntimes,timeindex,udiff)
      #})
      #nllk <- sum(nllk)
      #
     # break
    #}
    
    if(iteration<=1) likbase <- nllk
    if(iteration > maxit | 
       (iteration>2 & (nllk - likbase)/(oldlik - likbase) < 1 + tol) ) {
      nllk <- oldlik; break}
    
    #print(newdiff)
    if(print==TRUE){
        cat("iteration: ",iteration,#"; relative change: ", relchange,
            "; nllk: ",nllk,"\n")
        #print(l)
    }
    iteration <- iteration + 1
  }
  
  return(list(working_parm=parm,common_parm=parmbar[-indices],nllk=nllk))
}
