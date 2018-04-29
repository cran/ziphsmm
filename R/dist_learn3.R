
#######################################################
#' Distributed learning for a longitudinal continuous-time zero-inflated Poisson
#' hidden Markov model, where zero-inflation only happens in State 1 with covariates 
#' in the state-dependent parameters and transition rates.
#' @param ylist list of observed time series values for each subject
#' @param xlist list of design matrices for each subject. 
#' @param timelist list of time indices
#' @param M number of latent states
#' @param initparm matrix of initial working parameters for prior, transition, 
#' zero proportion, and emission parameters. 
#' @param yceil a scalar defining the ceiling of y, above which the values will be
#' truncated. Default to NULL. 
#' @param rho tuning parameter in the distributed learning algorithm. Default
#' to 1.
#' @param priorclust a vector to specify the grouping for state prior. Default to
#' NULL, which means no grouping.
#' @param tpmclust a vector to specify the grouping for the intercepts in state 
#' transition rates. Default to NULL, which means no grouping.
#' @param tpmslopeclust a vector to specify the grouping for the slopes in state 
#' transition rates. Default to NULL, which means no grouping.
#' @param emitclust a vector to specify the grouping for the intercepts in Poisson
#' regressions. Default to NULL, which means no grouping.
#' @param zeroclust a vector to specify the grouping for the intercepts in ZIP 
#' regression. Default to NULL, which means no grouping.
#' @param slopeclust a vector to specify the grouping for the slopes in Poisson and
#' ZIP regressions. Default to NULL, which means no grouping.
#' @param group a list containing group information.
#' @param maxit maximum number iteration. Default to 100.
#' @param tol tolerance in the terms of the relative change in the norm of the
#' common coefficients. Default to 1e-4. 
#' @param ncores number of cores to be used for parallel programming. Default to 1.
#' @param seed a seed for the random initialization of the algorithm
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
#' set.seed(12933)
#' nsubj <- 10
#' ns <- 2000
#' ylist <- vector(mode="list",length=nsubj)
#' xlist <- vector(mode="list",length=nsubj)
#' timelist <- vector(mode="list",length=nsubj)
#'
#' priorparm <- 0
#' tpmparm <- c(-2,0.1,-2,-0.2)
#' zeroindex <- c(1,0)
#' zeroparm <- c(0,0.5)
#' emitparm <- c(2,0.2,3,0.3)
#' workparm <- NULL
#' 
#' for(n in 1:nsubj){
#'  
#'  xlist[[n]] <- matrix(rep(c(0,1),rep(1000,2)),nrow=2000,ncol=1)
#'  
#'  timeindex <- rep(1,2000)
#'  for(i in 2:2000) timeindex[i] <- timeindex[i-1] + sample(1:4,1)
#'  timelist[[n]] <- timeindex
#'    
#'    workparm <- rbind(workparm,c(priorparm,tpmparm,zeroparm,emitparm))
#'  
#'  result <- hmmsim3.cont(workparm,2,2000,zeroindex,x=xlist[[n]],timeindex=timeindex)
#'  ylist[[n]] <- result$series
#' }
#'

#'
#' ####
#' M <- 2
#' priorclust <- c(rep(1,5),rep(2,5))
#' tpmclust <- c(rep(1,5),rep(2,5))
#' tpmslopeclust <- c(rep(1,5),rep(2,5))
#' zeroclust <- NULL
#' emitclust <- NULL
#' slopeclust <- rep(1,10)
#' 
#' group <- vector(mode="list",length=2)
#' group[[1]] <- 1:5; group[[2]] <- 6:10
#' ###

#' time <- proc.time()
#' result <- dist_learn3(ylist, xlist, timelist, 2,workparm, 
#'                      NULL, rho=0.1, priorclust,tpmclust,tpmslopeclust,
#'                      emitclust,zeroclust,slopeclust,group,ncores=1,
#'                      maxit=10, tol=1e-4, method="CG",print=TRUE)
#' proc.time() - time
#' }
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export

dist_learn3 <- function(ylist, xlist, timelist, M, initparm, yceil=NULL,
                        rho=1, priorclust=NULL,tpmclust=NULL,tpmslopeclust=NULL,
                        emitclust=NULL,zeroclust=NULL,slopeclust=NULL,group,
                        maxit=100, tol=1e-4, ncores=1,seed=0,
                        method="Nelder-Mead", print=TRUE, libpath=NULL,...){
  
  nsubj <- length(ylist)
  
  #retrieve working parameters
  ncolx <- ncol(xlist[[1]]) + 1
  allparm <- initparm
  
  #map [pi,tpmint,tpmslope,zeroint,zerocov,emit1int,emit1cov,...,emitMint,emitMcov]
  #to [pi,tpmint,tpmslope,zeroint,emit1int,...,emitMint,zerocov,emit1cov,...,emitMcov]
  mapf <- function(oldparm,M,ncolx){
    part1 <- oldparm[1:(M-1)]
    lid <- (M-1+M*(M-1)*ncolx)
    part2 <- oldparm[seq(M,length=M*(M-1),by=ncolx)]
    part3 <- oldparm[-c(1:(M-1),
                        seq(M,length=M*(M-1),by=ncolx),
                        (lid+1):length(oldparm))]
    
    part4 <- oldparm[seq(lid+1,length=M+1,by=ncolx)]
    part5 <- oldparm[-c(1:lid,seq(lid+1,length=M+1,by=ncolx))]
    return(c(part1,part2,part3,part4,part5))
  }
  
  invmapf <- function(newparm,M,ncolx){
    result <- rep(NA, M-1+M*(M-1)*ncolx+ncolx*(1+M))
    lid <- (M-1+M*(M-1)*ncolx)
    result[1:(M-1)] <- newparm[1:(M-1)]
    result[seq(M,length=M*(M-1),by=ncolx)] <- newparm[(M):(M+M*(M-1)-1)]
    result[-c(1:(M-1),
               seq(M,length=M*(M-1),by=ncolx),
               (lid+1):length(result))] <- newparm[(M-1+M*(M-1)+1):(M-1+M*(M-1)+M*(M-1)*(ncolx-1))]
    
    result[seq(lid+1,length=M+1,by=ncolx)] <- newparm[(lid+1):(lid+M+1)]
    
    result[-c(1:lid,seq(lid+1,length=M+1,by=ncolx))] <- newparm[(lid+M+2):length(result)]
    return(result)
  }
  
  #newparm <- mapf(allparm,M,ncolx)
  #invmapf(newparm,M,ncolx)
  ntotal <- ncol(allparm)
  
  #########################
  #initial J matrix
  J <- diag(1,ntotal)
  lz <- 0
  #cannot be totally subject-specific
  totalgroup <- max(c(priorclust,tpmclust,tpmslopeclust,emitclust,zeroclust,slopeclust,1))
  commonindex <- NULL #for z
  last <- 0
  lastcommon <- 0
  lastcluster <- 0
  clusterindex <- vector(mode="list",length=totalgroup)
  
  last2 <- 0 #for theta
  lastcommon2 <- 0
  lastcluster2 <- 0
  commonindex2 <- NULL
  clusterindex2 <- NULL
  
  jcommonindex <- NULL #for constraint id in J
  jclusterindex <- NULL
  lastjcommon <- 0
  lastjcluster <- 0
  lastj <- 0
  rowtodelete <- NULL
  
  if(is.null(priorclust)){#subject specific
    rowtodelete <- c(rowtodelete,1:(M-1))}else if(max(priorclust)==1){#common
      lz <- lz + (M-1)*max(priorclust)
      commonindex <- c(commonindex,seq(1,length=M-1,by=1))
      lastcommon <- max(commonindex)
      commonindex2 <- c(commonindex2,1:(M-1))
      lastcommon2 <- M-1
      jcommonindex <- c(jcommonindex,1:(M-1))
      lastjcluster <- M-1
    }else{#clustering
      lz <- lz + (M-1)*max(priorclust)
      for(g in 1:totalgroup) clusterindex[[g]] <- c(clusterindex[[g]],
                                                    seq(g,length=M-1,by=totalgroup))
      lastcluster <- max(clusterindex[[totalgroup]])
      clusterindex2 <- c(clusterindex2,1:(M-1))
      lastcluster2 <- M-1
      jclusterindex <- c(jclusterindex,1:(M-1))
      lastjcluster <- M-1
    }
  last <- max(c(lastcommon,lastcluster))
  last2 <- max(c(lastcommon2,lastcluster2))
  lastj <- max(lastjcluster,lastjcommon)
  
  if(is.null(tpmclust)){#subject specific
    rowtodelete <- c(rowtodelete,M:(M*M-1))}else if(max(tpmclust)==1){#common
      lz <- lz + M*(M-1)*max(tpmclust)
      commonindex <- c(commonindex,seq(last+1,length=M*(M-1),by=1))
      lastcommon <- max(commonindex)
      commonindex2 <- c(commonindex2,M:(M*M-1))
      lastcommon2 <- M*M-1
      jcommonindex <- c(jcommonindex,(lastj+1):(lastj+M*(M-1)))
      lastjcommon <- max(jcommonindex)
    }else{#clustering
      lz <- lz + M*(M-1)*max(tpmclust)
      for(g in 1:totalgroup) clusterindex[[g]] <- c(clusterindex[[g]],
                                                    seq(last+g,length=M*(M-1),by=totalgroup))
      lastcluster <- max(clusterindex[[totalgroup]])
      clusterindex2 <- c(clusterindex2,M:(M*M-1))
      lastcluster2 <- M*M-1
      jclusterindex <- c(jclusterindex,(lastj+1):(lastj+M*(M-1)))
      lastjcluster <- max(jclusterindex)
    }
  last <- max(c(lastcommon,lastcluster))
  last2 <- max(c(lastcommon2,lastcluster2))
  lastj <- max(lastjcluster,lastjcommon)
  
  
  if(is.null(tpmslopeclust)){#subject specific
    rowtodelete <- c(rowtodelete,(M*M):(M*M+M*(M-1)*(ncolx-1)-1))}else if(max(tpmslopeclust)==1){#common
      lz <- lz + M*(M-1)*max(tpmslopeclust)
      commonindex <- c(commonindex,seq(last+1,length=M*(M-1)*(ncolx-1),by=1))
      lastcommon <- max(commonindex)
      commonindex2 <- c(commonindex2,(M*M):(M*M+M*(M-1)*(ncolx-1)-1))
      lastcommon2 <- M*M+M*(M-1)*(ncolx-1)-1
      jcommonindex <- c(jcommonindex,(lastj+1):(lastj+M*(M-1)*(ncolx-1)))
      lastjcommon <- max(jcommonindex)
    }else{#clustering
      lz <- lz + M*(M-1)*(ncolx-1)*max(tpmslopeclust)
      for(g in 1:totalgroup) clusterindex[[g]] <- c(clusterindex[[g]],
                                                    seq(last+g,length=M*(M-1),by=totalgroup))
      lastcluster <- max(clusterindex[[totalgroup]])
      clusterindex2 <- c(clusterindex2,(M*M):(M*M+M*(M-1)*(ncolx-1)-1))
      lastcluster2 <- M*M+M*(M-1)*(ncolx-1)-1
      jclusterindex <- c(jclusterindex,(lastj+1):(lastj+M*(M-1)*(ncolx-1)))
      lastjcluster <- max(jclusterindex)
    }
  last <- max(c(lastcommon,lastcluster))
  last2 <- max(c(lastcommon2,lastcluster2))
  lastj <- max(lastjcluster,lastjcommon)
  
  
  if(is.null(zeroclust)){#subject specific
    rowtodelete <- c(rowtodelete,(M*M+M*(M-1)*(ncolx-1))) }else if(max(zeroclust)==1){#common
      lz <- lz + max(zeroclust)
      commonindex <- c(commonindex,last+1)
      lastcommon <- max(commonindex)
      commonindex2 <- c(commonindex2,(M*M+M*(M-1)*(ncolx-1)))
      lastcommon2 <- (M*M+M*(M-1)*(ncolx-1))
      jcommonindex <- c(jcommonindex,lastj+1)
      lastjcommon <- lastj+1
    }else{#clustering
      lz <- lz + max(zeroclust)
      for(g in 1:totalgroup) clusterindex[[g]] <- c(clusterindex[[g]],last+g)
      lastcluster <- max(clusterindex[[totalgroup]])
      clusterindex2 <- c(clusterindex2,(M*M+M*(M-1)*(ncolx-1)))
      lastcluster2 <- (M*M+M*(M-1)*(ncolx-1))
      jclusterindex <- c(jclusterindex,lastj+1)
      lastjcluster <- lastj+1
    }
  last <- max(c(lastcommon,lastcluster))
  last2 <- max(c(lastcommon2,lastcluster2))
  lastj <- max(c(lastjcommon,lastjcluster))
  
  if(is.null(emitclust)){#subject specific
    rowtodelete <- c(rowtodelete,(M*M+M*(M-1)*(ncolx-1)+1):(M*M+M*(M-1)*(ncolx-1)+M)) }else if(max(emitclust)==1){#common
      lz <- lz + M*max(emitclust)
      commonindex <- c(commonindex,seq(last+1,length=M,by=1))
      lastcommon <- max(commonindex)
      commonindex2 <- c(commonindex2,(M*M+M*(M-1)*(ncolx-1)+1):(M*M+M*(M-1)*(ncolx-1)+M))
      lastcommon2 <- max(commonindex2)
      jcommonindex <- c(jcommonindex,(lastj+1):(lastj+M))
      lastjcommon <- max(jcommonindex)
    }else{
      lz <- lz + M*max(emitclust)
      for(g in 1:totalgroup) clusterindex[[g]] <- c(clusterindex[[g]],
                                                    seq(last+g,length=M,by=totalgroup))
      lastcluster <- max(clusterindex[[totalgroup]])
      clusterindex2 <- c(clusterindex2,(M*M+M*(M-1)*(ncolx-1)+1):(M*M+M*(M-1)*(ncolx-1)+M))
      lastcluster2 <- max(clusterindex2)
      jclusterindex <- c(jclusterindex,(lastj+1):(lastj+M))
      lastjcluster <- max(jclusterindex)
    }
  last <- max(c(lastcommon,lastcluster))
  last2 <- max(c(lastcommon2,lastcluster2))
  lastj <- max(c(lastjcommon,lastjcluster))
  
  if(is.null(slopeclust)){#subject specific
    rowtodelete <- c(rowtodelete,(M*M+M*(M-1)*(ncolx-1)+M+1):ntotal) }else if(max(slopeclust)==1){#common
      lz <- lz + (ncolx-1)*(M+1)*max(slopeclust)
      commonindex <- c(commonindex,seq(last+1,length=(ncolx-1)*(M+1),by=1))
      lastcommon <- max(commonindex)
      commonindex2 <- c(commonindex2,(M*M+M*(M-1)*(ncolx-1)+M+1):ntotal)
      lastcommon2 <- max(commonindex2)
      jcommonindex <- c(jcommonindex,(lastj+1):(lastj+(ncolx-1)*(M+1)))
      lastjcommon <- max(jcommonindex)
    }else{
      lz <- lz + (ncolx-1)*(M+1)*max(slopeclust)
      for(g in 1:totalgroup) clusterindex[[g]] <- c(clusterindex[[g]],
                                                    seq(last+g,length=(ncolx-1)*(M+1),by=totalgroup))
      lastcluster <- max(clusterindex[[totalgroup]])
      clusterindex2 <- c(clusterindex2,(M*M+M*(M-1)*(ncolx-1)+M+1):ntotal)
      lastcluster2 <- max(clusterindex2)
      jclusterindex <- c(jclusterindex,(lastj+1):(lastj+(ncolx-1)*(M+1)))
      lastjcluster <- max(jclusterindex)
    }
  last <- max(c(lastcommon,lastcluster))
  last2 <- max(c(lastcommon2,lastcluster2))
  lastj <- max(c(lastjcommon,lastjcluster))
  
  
  #paramters and their gradients
  if(!is.null(rowtodelete)) J <- J[-rowtodelete,]
  
  set.seed(seed)
  parm <- t(sapply(1:nrow(allparm), function(gg){
    thisrand <- runif(ncol(allparm),-0.01,0.01)
    this <- allparm[gg,] + thisrand
    mapf(this,M,ncolx)
  }))
  
  l <- matrix(0,nsubj,lz)
  z <- numeric(lz)
  
  
  #must have some common effects
  #otherwise, just split into subgroups and refit
  #must have some common effects
  #otherwise, just split into subgroups and refit
  if(is.null(commonindex)){
    print("Must have some common effects! Otherwise, simply split into clusters and refit.")
  }else if(totalgroup>1){ #some clustering some common
    #most common case
    
    tempcluster <- vector(mode="list",length=totalgroup)
    olddiff <- 0
    #initial value
    for(g in 1:totalgroup) {
      tempcluster[[g]] <- 
        J[jclusterindex,clusterindex2]%*%t(parm[,clusterindex2])[,group[[g]]]
      z[clusterindex[[g]]] <- rowMeans(tempcluster[[g]]) + 
        colMeans(l[group[[g]],clusterindex[[g]]])/rho
      #primal residual
      olddiff <- olddiff + sum((tempcluster[[g]]-z[clusterindex[[g]]])^2)
    }
    tempcommon <- J[jcommonindex,commonindex2]%*%t(parm[,commonindex2])
    z[commonindex] <- rowMeans(tempcommon) + colMeans(l[,commonindex])/rho
    #primal residual
    olddiff <- olddiff + sum((tempcommon-z[commonindex])^2)
    oldnorm <- sum(z^2)
    
    olddualdiff <- 0
    olddualparm <- l
    olddualnorm <- 0
    
    #new functions for penalized negloglik
    #new functions for penalized negloglik
    zipnegloglik_cov_cont3 <- ziphsmm::zipnegloglik_cov_cont3
    newf <- function(initparm,y,covariates,M,ntimes,timeindex,
                     zi,rho,li){
      cov <- cbind(1,covariates)
      part1 <- zipnegloglik_cov_cont3(initparm,y,cov,M,ntimes,timeindex)
      parmnew <- mapf(initparm,M,ncol(covariates)+1)
      diff <- J%*%parmnew - zi
      part2 <- t(li)%*%diff
      part3 <- 0.5*rho*t(diff)%*%diff
      return(part1+part2+part3)
    }
    
    
    newgradf <- function(initparm,y,covariates,M,ntimes,timeindex,
                         zi,rho,li){
      cov <- cbind(1,covariates)
      part1 <- pracma::grad(zipnegloglik_cov_cont3,initparm,
                              y=y,covariates=covariates,M=M,
                              ntimes=ntimes,timeindex=timeindex)
      part2 <- t(J)%*%li
      parmnew <- mapf(initparm,M,ncol(covariates)+1)
      part3 <- rho * (t(J) %*% ( J%*%parmnew - zi))
      return(part1+part2+part3)
    }
    
    #start iterations
    iteration <- 1
    nllk <- 0
    dual_change <- NULL
    nllk_change <- NULL
    primal_change <- NULL
    resid_change <- NULL
    
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
          x <- xlist[[i]]
          if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
          timeindex <- timelist[[i]]
          ntimes <- length(y)
        
          #get subject-specific z and l
          for(kk in 1:totalgroup)
            if(i%in%group[[kk]]){gi <- kk}else{next}
          fullindex <-  c(clusterindex[[gi]],commonindex)  #for z and l
          zi <- z[fullindex]
          li <- l[i,fullindex]
          initparm <- invmapf(parm[i,],M,ncolx)
          
          optim(par=initparm,fn=newf,gr=newgradf,
                M=M,y=y,covariates=x,ntimes=ntimes,
                timeindex=timeindex,
                zi=zi,rho=newrho,li=li,
                method=method,...)
          #newf(initparm,y,x,M,ntimes,timeindex,zi,rho,li)
          #newgradf(initparm,y,x,M,ntimes,timeindex,zi,rho,li)
        })
        # proc.time() - time
      }else{
        
        cl <- parallel::makeCluster(ncores)
        parallel::clusterExport(cl,c("M","ylist","xlist","timelist","yceil","l","parm",
                                     "z","newrho","method","mapf","invmapf","group",
                                     "newf","newgradf","grad_zipnegloglik_nocov_cont",
                                     "zipnegloglik_nocov_cont","J","libpath",
                                     "totalgroup","clusterindex","commonindex"),
                                envir=environment())
        #time <- proc.time()
        tempresult <- parallel::parLapply(cl, 1:nsubj, function(i){
          if(!is.null(libpath)) .libPaths(libpath)  #'~/R_p4/library'
          library(ziphsmm)
          y <- ylist[[i]]
          x <- xlist[[i]]
          if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
          timeindex <- timelist[[i]]
          ntimes <- length(y)
          
 
          for(kk in 1:totalgroup)
            if(i%in%group[[kk]]){gi <- kk}else{next}
          fullindex <-  c(clusterindex[[gi]],commonindex)  #for z and l
          zi <- z[fullindex]
          li <- l[i,fullindex]
          initparm <- invmapf(parm[i,],M,ncolx)
          
          optim(par=initparm,fn=newf,gr=newgradf,
                M=M,y=y,covariates=x,ntimes=ntimes,
                timeindex=timeindex, 
                zi=zi,rho=newrho,li=li,
                method=method)
          
        })
        parallel::stopCluster(cl)
        #proc.time()-time
      }
      
      #############################################################
      nllk <- sum(sapply(1:nsubj,function(i)tempresult[[i]]$value))
      parm <- t(sapply(1:nsubj,function(i) {
        temppar <- tempresult[[i]]$par
        
        mapf(temppar,M,ncolx)
      }))
      
      #update z
      newdiff <- 0
      for(g in 1:totalgroup) {
        tempcluster[[g]] <- 
          J[jclusterindex,clusterindex2]%*%t(parm[,clusterindex2])[,group[[g]]]
        z[clusterindex[[g]]] <- rowMeans(tempcluster[[g]]) + 
          colMeans(l[group[[g]],clusterindex[[g]]])/rho
        #primal residual
        newdiff <- newdiff + sum((tempcluster[[g]]-z[clusterindex[[g]]])^2)
      }
      tempcommon <- J[jcommonindex,commonindex2]%*%t(parm[,commonindex2])
      z[commonindex] <- rowMeans(tempcommon) + colMeans(l[,commonindex])/rho
      newdiff <- newdiff + sum((tempcommon-z[commonindex])^2)
      
      resid <- abs(newdiff-olddiff) / (1+olddiff)  
      resid_change <- c(resid_change, resid)
      
      newnorm <- sum(z^2)
      primal_diff <- abs(sqrt(newnorm) - sqrt(oldnorm)) / (1+sqrt(oldnorm))
      primal_change <- c(primal_change, primal_diff)
       
      
      #update l
      for(i in 1:nsubj){
        for(kk in 1:totalgroup)
          if(i%in%group[[kk]]){gi <- kk}else{next}
        fullindex <-  c(clusterindex[[gi]],commonindex)  #for z and l
        l[i,fullindex] <- l[i,fullindex] + rho * 
          (J%*%parm[i,]-z[fullindex])
      }
      
      newdualparm <- l
      newdualnorm <- sum(l^2)
      newdualdiff <- sum((newdualparm - olddualparm)^2)
      reldualchange <- sqrt(newdualdiff) / (sqrt(olddualnorm) +1)
      dual_change <- c(dual_change, reldualchange)
      
      if(iteration<=1) likbase <- nllk
      new_nllk_change <- abs(nllk-oldlik)/(1+oldlik)
      nllk_change <- c(nllk_change,new_nllk_change)
      
      kkt_cur <- max(primal_diff, reldualchange, new_nllk_change)#newzdiff
      if(iteration > maxit | 
         (iteration>2 & kkt_cur < tol )) {
        nllk <- oldlik; break}
      
      if(print==TRUE & iteration>=2){
        cat("iter:",iteration, "; threshold:", kkt_cur,"\n")
      }
      
      olddiff <- newdiff #
      
      olddualdiff <- newdualdiff
      olddualparm <- newdualparm
      olddualnorm <- newdualnorm
      oldnorm <- newnorm
      old_nllk_change <- new_nllk_change
      iteration <- iteration + 1
    }
    #reorder back
    workingparm <- t(sapply(1:nrow(parm),function(kkk) invmapf(parm[kkk,],M,ncolx)))
    return(list(working_parm=workingparm,
    #common_parm=z,
                threshold=pmax(primal_change[-1],dual_change[-1],nllk_change[-1]),
                nllk=nllk))
  }
  
}


