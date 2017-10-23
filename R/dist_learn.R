
#######################################################
#' Distributed learning for a longitudinal continuous-time zero-inflated Poisson
#' hidden Markov model, where zero-inflation only happens in State 1. Assume that
#' prior, transition rate are subject-specific, but the state-dependent parameters
#' are the same across subjects. A distributed learning algorithm
#' is used to estimate the parameters.
#' @param ylist list of observed time series values for each subject
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
#' @param ncores number of cores to be used for parallel programming. Default to 1.
#' @param method method for the distributed optimization in the ADMM framework.
#' @param print whether to print each iteration. Default to TRUE.
#' @param libpath path for the ziphsmm library if not the default set up. Default to NULL.
#' @param ... Further arguments passed on to the optimization methods
#' @return the maximum likelihood estimates of the zero-inflated hidden Markov model
#' @references Boyd, S., Parikh, N., Chu, E., Peleato, B. and Eckstein, J., 2011. 
#' Distributed optimization and statistical learning via the alternating direction method 
#' of multipliers. Foundations and Trends in Machine Learning, 3(1), pp.1-122.
#' @examples
#' \dontrun{
#' set.seed(930518)
#' nsubj <- 20
#' ns <- 10800
#' ylist <- vector(mode="list",length=nsubj)
#' timelist <- vector(mode="list",length=nsubj)

#' for(n in 1:nsubj){
#'   rn1 <- runif(1,-0.02,0.02)
#'   rn2 <- runif(1,-0.02,0.02)
#'  prior <- c(0.5 + rn1,0.2 + rn2,0.3 - rn1 - rn2)
#'  omega <- matrix(c(-0.3-rn1-rn2,0.2+rn1,0.1+rn2,
#'                    0.1+rn1,-0.2-rn1-rn2,0.1+rn2,
#'                    0.15+rn1,0.2+rn2,-0.35-rn1-rn2),3,3,byrow=TRUE)
#'  emit <- c(50,200,600)
#'  zero <- c(0.2,0,0)
#'  timeindex <- rep(1,1440)
#'  for(i in 2:ns) timeindex[i] <- timeindex[i-1] + sample(1:4,1)
#'  timelist[[n]] <- timeindex
#'  
#'  result <- hmmsim.cont(ns, 3, prior, omega, emit, zero, timeindex)
#'  ylist[[n]] <- result$series
#' }
#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(50, 225, 650)
#' zero_init <- 0.2
#' tpm_init <- matrix(c(-0.3,0.2,0.1,0.1,-0.2,0.1,0.15,0.2,-0.35),3,3,byrow=TRUE)
#'
#' result <- dist_learn(ylist, timelist, prior_init, tpm_init, 
#'                     emit_init, zero_init, rho=10, 
#'                     maxit=100, tol=1e-4, method="CG")
#' }
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export

dist_learn <- function(ylist, timelist, prior_init, tpm_init, 
                        emit_init, zero_init, yceil=NULL,
                        rho=1, maxit=100, tol=1e-4, ncores=1,
                        method="CG", print=TRUE, libpath=NULL, ...){
  
  nsubj <- length(ylist)
  M <- ncol(tpm_init)
  
  #retrieve working parameters
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
  
  
  #fixed indices in working parameters
  indices <- 1:(M*M-1)
  #others
  ntotal <- length(allparm)
  nspecial <- length(indices)
  common_indices <- setdiff(1:ntotal, indices)
  
  J <- matrix(0, nrow=ntotal-nspecial, ncol=ntotal)
  for(i in 1:(ntotal-nspecial)) J[i,common_indices[i]] <- 1
  
  
  #paramters and their gradients
  parm <- matrix(0, nsubj, length(allparm))
  #set.seed(518930)
  for(i in 1:nsubj) parm[i,] <- allparm #+ c(rep(0,nspecial),
                                           # runif(ntotal-nspecial,-0.05,0.05))
  
  parmbar <- colMeans(parm)
  olddiff <- parm[,-indices,drop=FALSE]-
    matrix(rep(parmbar[-indices],nsubj),
           nsubj,length(allparm)-length(indices),byrow=TRUE)
  olddiff <- sum(olddiff^2)
  oldnorm <- sum(parm[,-indices]^2)
  olddiff <- olddiff / (1+oldnorm)
  #olddiff <- parm-matrix(rep(parmbar,nsubj),nrow=nsubj,byrow=TRUE)
  #oldparm <- c(as.vector(parm),parmbar[-indices]) #
  #oldnorm <- sum(oldparm^2) #
  
  grad <- matrix(0, nsubj, length(allparm))
  l <- matrix(0,nsubj,ntotal - nspecial)
  
  olddualdiff <- 0
  olddualparm <- l
  olddualnorm <- sum(l^2)
  #ADMM
  #thetai := argmin(fi(thetai) + t(li)%*%(thetai-thetabar) + 0.5*rho*|thetai-thetabar|^2)
  #thetabar := mean(thetai)
  #li := li + rho(thetai - thetabar)
  zipnegloglik_nocov_cont <- ziphsmm::zipnegloglik_nocov_cont
  newf <- function(initparm,y,M,ntimes,timeindex,udiff,
                   parmbar,indices,rho,li){
    part1 <- zipnegloglik_nocov_cont(initparm,M, y,ntimes,timeindex,udiff)
    diff <- initparm[-indices] - parmbar[-indices]
    part2 <- t(li)%*%diff
    part3 <- 0.5*rho*t(diff)%*%diff
    return(part1+part2+part3)
  }
  
  grad_zipnegloglik_nocov_cont <- ziphsmm::grad_zipnegloglik_nocov_cont
  newgradf <- function(initparm,y,M,ntimes,timeindex,udiff,
                       parmbar,indices,rho,li){
    part1 <- grad_zipnegloglik_nocov_cont(initparm,M, y,ntimes,timeindex,udiff)
    part2 <- t(J)%*%li
    part3 <- rho * (t(J)%*%J) %*% (initparm - parmbar)
    return(part1+part2+part3)
  }
  
  iteration <- 1
  nllk <- 0
  primal_change <- NULL
  dual_change <- NULL
  nllk_change <- NULL
  
  
  #recursion
  while(iteration<=maxit){
    newrho <- rho
    #newrho <- rho * iteration^(-1)
    #distributed
    oldlik <- nllk
    if(ncores==1){
      #time <- proc.time()
        tempresult <- lapply(1:nsubj, function(i){
          y <- ylist[[i]]
          if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
          timeindex <- timelist[[i]]
          ntimes <- length(y)
          
          vdiff <- diff(timeindex)
          udiff <- sort(unique(vdiff))
          li <- l[i,]
          initparm <- parm[i,]
          
          optim(par=initparm,fn=newf,gr=newgradf,
                M=M,y=y,ntimes=ntimes,timeindex=timeindex,udiff=udiff, 
                parmbar=parmbar,indices=indices,rho=newrho,li=li,
                method=method,...)
          #newf(initparm,y,M,ntimes,timeindex,udiff,parmbar,indices,rho,li)
          #grad(newf,initparm,M=M,y=y,ntimes=ntimes,timeindex=timeindex,udiff=udiff, 
          #   parmbar=parmbar,indices=indices,rho=newrho,li=li)
          #newgradf(initparm,y,M,ntimes,timeindex,udiff,parmbar,indices,rho,li)
        })
       # proc.time() - time
    }else{
     
      cl <- parallel::makeCluster(ncores)
      parallel::clusterExport(cl,c("M","ylist","timelist","yceil","l","parm",
                         "indices","parmbar","newrho","method",
                         "newf","newgradf","grad_zipnegloglik_nocov_cont",
                         "zipnegloglik_nocov_cont","J","libpath"),envir=environment())
      #time <- proc.time()
      tempresult <- parallel::parLapply(cl, 1:nsubj, function(i){
        if(!is.null(libpath)) .libPaths(libpath)  #'~/R_p4/library'
        library(ziphsmm)
        y <- ylist[[i]]
        if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
        timeindex <- timelist[[i]]
        ntimes <- length(y)
        
        vdiff <- diff(timeindex)
        udiff <- sort(unique(vdiff))
        li <- l[i,]
        initparm <- parm[i,]
        
        optim(par=initparm,fn=newf,gr=newgradf,
              M=M,y=y,ntimes=ntimes,timeindex=timeindex,udiff=udiff, 
              parmbar=parmbar,indices=indices,rho=newrho,li=li,
              method=method)
         
      })
      parallel::stopCluster(cl)
      #proc.time()-time
    }
    
    #############################################################
    nllk <- sum(sapply(1:nsubj,function(i)tempresult[[i]]$value))
    #permutation of states
    parm <- t(sapply(1:nsubj,function(i) {
      temppar <- tempresult[[i]]$par
      c(temppar[1:(M*M)], sort(temppar[(M*M+1):length(temppar)]))
      }))
    
    parmbar <- colMeans(parm)
    
    for(i in 1:nsubj)
      l[i,] <- l[i,] + newrho*(parm[i,-indices] - parmbar[-indices])
    
    newdiff <- parm[,-indices,drop=FALSE]-
      matrix(rep(parmbar[-indices],nsubj),
             nsubj,ntotal-nspecial,byrow=TRUE)#+
     # matrix(rep(colMeans(l),nsubj),
     #        nsubj,ntotal-nspecial,byrow=TRUE)/newrho
    newdiff <- sum(newdiff^2)
    newnorm <- sum(parm[,-indices]^2)
    relchange <- newdiff / (1+oldnorm)
    #newparm <- c(as.vector(parm),parmbar[-indices])
    #newnorm <- sum(newparm^2)
    #newdiff <- sum((newparm - oldparm)^2)
    #relchange <- newdiff / (oldnorm +1)
    primal_change <- c(primal_change,relchange)
    
    newdualparm <- l
    newdualnorm <- sum(l^2)
    newdualdiff <- sum((newdualparm - olddualparm)^2)
    reldualchange <- newdualdiff / (olddualnorm +1)
    dual_change <- c(dual_change, reldualchange)
    
    if(iteration<=1) likbase <- nllk
    new_nllk_change <- abs(nllk-oldlik)/(1+oldlik)
    nllk_change <- c(nllk_change,new_nllk_change)
    
    if(iteration > maxit | 
       (iteration>2 & relchange < tol & reldualchange < tol )) {
      nllk <- oldlik; break}
    
    if(print==TRUE & iteration>=2){
      cat("iter:",iteration, "; residual_rel_chg:", relchange,
          "; dual_rel_chg:",reldualchange,"; nllk_rel_chg:",new_nllk_change,"\n")
    }
    
    olddiff <- newdiff #
    #oldparm <- newparm #
    oldnorm <- newnorm #
    olddualdiff <- newdualdiff
    olddualparm <- newdualparm
    olddualnorm <- newdualnorm
    
    old_nllk_change <- new_nllk_change
    iteration <- iteration + 1
  }
  
  #natparm <- vector(mode="list",length=nsubj)
  #natparm <- retrieve_nocov_cont(fitparm,M)
  
  return(list(working_parm=parm,common_parm=parmbar[-indices],
              residual_change=primal_change[-1],dual_change=dual_change[-1],
              nllk=nllk))
}



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

#if(iteration<=1) likbase <- nllk
#if(iteration > maxit | 
#   (iteration>2 & (nllk - likbase)/(oldlik - likbase) < 1 + tol) ) {
#  nllk <- oldlik; break}

