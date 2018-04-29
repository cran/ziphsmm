
#######################################################
#' Distributed learning for a longitudinal continuous-time zero-inflated Poisson
#' hidden Markov model, where zero-inflation only happens in State 1. Assume that
#' priors, transition rates and state-dependent parameters can be subject-specific,
#' clustered by group, or common. But at least one set of the parameters have to be 
#' common across all subjects.
#' @param ylist list of observed time series values for each subject
#' @param timelist list of time indices
#' @param prior_init a vector of initial values for prior probability for each state
#' @param tpm_init a matrix of initial values for transition rate matrix
#' @param emit_init a vector of initial values for the means for each poisson distribution
#' @param zero_init a scalar initial value for the structural zero proportion 
#' @param yceil a scalar defining the ceiling of y, above which the values will be
#' truncated. Default to NULL. 
#' @param rho tuning parameters in the distributed learning algorithm. Default to 1.
#' @param priorclust a vector to specify the grouping for state prior. Default to
#' NULL, which means no grouping.
#' @param tpmclust a vector to specify the grouping for state transition rates. 
#' Default to NULL, which means no grouping.
#' @param emitclust a vector to specify the grouping for Poisson means. Default to
#' NULL, which means no grouping.
#' @param zeroclust a vector to specify the grouping for structural zero proportions. 
#' Default to NULL, which means no grouping.
#' @param group a list containing group information.
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
#' nsubj <- 10
#' ns <- 5040
#' ylist <- vector(mode="list",length=nsubj)
#' timelist <- vector(mode="list",length=nsubj)

#' prior1 <- c(0.5,0.2 ,0.3 )
#' omega1 <- matrix(c(-0.3,0.2,0.1,
#'                   0.1,-0.2,0.1,
#'                   0.15,0.2,-0.35),3,3,byrow=TRUE)
#' prior2 <- c(0.3,0.3 ,0.4 )
#' omega2 <- matrix(c(-0.5,0.25,0.25,
#'                    0.2,-0.4,0.2,
#'                    0.15,0.3,-0.45),3,3,byrow=TRUE)
#' emit <- c(50,200,600)
#' zero <- c(0.2,0,0)

#' for(n in 1:nsubj){
  
#'  timeindex <- rep(1,ns)
#'  for(i in 2:ns) timeindex[i] <- timeindex[i-1] + sample(1:4,1)
#'  timelist[[n]] <- timeindex
  
#'  if(n<=5){
#'    result <- hmmsim.cont(ns, 3, prior1, omega1, emit, zero, timeindex)
#'    ylist[[n]] <- result$series
#'  }else{
#'    result <- hmmsim.cont(ns, 3, prior2, omega2, emit, zero, timeindex)
#'    ylist[[n]] <- result$series
#'  }
#' }

#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(50, 225, 650)
#' zero_init <- 0.2
#' tpm_init <- matrix(c(-0.3,0.2,0.1,0.1,-0.2,0.1,0.15,0.2,-0.35),3,3,byrow=TRUE)
 
#' M <- 3
#' priorclust <- NULL
#' tpmclust <- c(1,1,1,1,1,2,2,2,2,2)
#' zeroclust <- rep(1,10)
#' emitclust <- rep(1,10)
#' group <- vector(mode="list",length=2)
#' group[[1]] <- 1:5; group[[2]] <- 6:10
#' result <- dist_learn(ylist, timelist, prior_init, tpm_init, 
#'                     emit_init, zero_init,NULL, rho=1,priorclust,tpmclust,
#'                     emitclust,zeroclust,group,ncores=1,
#'                     maxit=50, tol=1e-4, method="CG", print=TRUE)
#' }
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export

dist_learn <- function(ylist, timelist, prior_init, tpm_init, 
                        emit_init, zero_init, yceil=NULL,
                        rho=1, priorclust=NULL,tpmclust=NULL,
                        emitclust=NULL,zeroclust=NULL,group,
                        maxit=100, tol=1e-4, ncores=1,
                        method="Nelder-Mead", print=TRUE, libpath=NULL, ...){
  
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
  ntotal <- length(allparm)
  
  #########################
  #initial J matrix
  J <- diag(1,ntotal)
  lz <- 0
  #cannot be totally subject-specific
  
  totalgroup <- max(c(priorclust,tpmclust,emitclust,zeroclust,1))
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
  
  rowtodelete <- NULL  #row to delete in j
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
  
  if(is.null(zeroclust)){#subject specific
    rowtodelete <- c(rowtodelete,M*M) }else if(max(zeroclust)==1){#common
      lz <- lz + max(zeroclust)
      commonindex <- c(commonindex,last+1)
      lastcommon <- max(commonindex)
      commonindex2 <- c(commonindex2,M*M)
      lastcommon2 <- M*M
      jcommonindex <- c(jcommonindex,lastj+1)
      lastjcommon <- lastj+1
    }else{#clustering
      lz <- lz + max(zeroclust)
      for(g in 1:totalgroup) clusterindex[[g]] <- c(clusterindex[[g]],last+g)
      lastcluster <- max(clusterindex[[totalgroup]])
      clusterindex2 <- c(clusterindex2,M*M)
      lastcluster2 <- M*M
      jclusterindex <- c(jclusterindex,lastj+1)
      lastjcluster <- lastj+1
    }
  last <- max(c(lastcommon,lastcluster))
  last2 <- max(c(lastcommon2,lastcluster2))
  lastj <- max(c(lastjcommon,lastjcluster))
  
  if(is.null(emitclust)){#subject specific
    rowtodelete <- c(rowtodelete,(M*M+1):(M*M+M)) }else if(max(emitclust)==1){#common
      lz <- lz + M*max(emitclust)
      commonindex <- c(commonindex,seq(last+1,length=M,by=1))
      commonindex2 <- c(commonindex2,(M*M+1):(M*M+M))
      jcommonindex <- c(jcommonindex,(lastj+1):(lastj+M))
    }else{
      lz <- lz + M*max(emitclust)
      for(g in 1:totalgroup) clusterindex[[g]] <- c(clusterindex[[g]],
                                            seq(last+g,length=M,by=totalgroup))
      clusterindex2 <- c(clusterindex2,(M*M+1):(M*M+M))
      jclusterindex <- c(jclusterindex,(lastj+1):(lastj+M))
    }
  
  #paramters and their gradients
  if(!is.null(rowtodelete)) J <- J[-rowtodelete,]
  parm <- matrix(0, nsubj, length(allparm))
  #set.seed(518930)
  for(i in 1:nsubj) parm[i,] <- allparm + runif(ntotal,-0.05,0.05)
  
  l <- matrix(0,nsubj,lz)
  z <- numeric(lz)
  
  
  
  #must have some common effects
  #otherwise, just split into subgroups and refit
  if(is.null(commonindex)) {
    print("Must have some common effects! Otherwise, simply split into clusters and refit.")
  }else if(totalgroup==1){#some common, some subject-specific
    olddiff <- 0
    tempcommon <- J[jcommonindex,commonindex2]%*%t(parm[,commonindex2])
    z[commonindex] <- rowMeans(tempcommon) + colMeans(l[,commonindex])/rho
    
    #primal residual
    olddiff <- sum((tempcommon-z[commonindex])^2)
    #
    oldnorm <- sum(z^2)
    
    olddualdiff <- 0
    olddualparm <- l
    olddualnorm <- 0
    
    #new functions for penalized negloglik
    zipnegloglik_nocov_cont <- ziphsmm::zipnegloglik_nocov_cont
    newf <- function(initparm,y,M,ntimes,timeindex,udiff,
                     zi,rho,li){
      part1 <- zipnegloglik_nocov_cont(initparm,M, y,ntimes,timeindex,udiff)
      diff <- J%*%initparm - zi
      part2 <- t(li)%*%diff
      part3 <- 0.5*rho*t(diff)%*%diff
      return(part1+part2+part3)
    }
    
    grad_zipnegloglik_nocov_cont <- ziphsmm::grad_zipnegloglik_nocov_cont
    newgradf <- function(initparm,y,M,ntimes,timeindex,udiff,
                         zi,rho,li){
      part1 <- grad_zipnegloglik_nocov_cont(initparm,M, y,ntimes,timeindex,udiff)
      part2 <- t(J)%*%li
      part3 <- rho * (t(J) %*% ( J%*%initparm - zi))
      return(part1+part2+part3)
    }
    
    #start iterations
    iteration <- 1
    nllk <- 0
    resid_change <- NULL
    dual_change <- NULL
    nllk_change <- NULL
    primal_change <- NULL
    
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
          
          fullindex <-  commonindex  #for z and l
          zi <- z[fullindex]
          li <- l[i,fullindex]
          initparm <- parm[i,]
          
          optim(par=initparm,fn=newf,gr=newgradf,
                M=M,y=y,ntimes=ntimes,timeindex=timeindex,udiff=udiff, 
                zi=zi,rho=newrho,li=li,
                method=method,...)
          #newf(initparm,y,M,ntimes,timeindex,udiff,zi,rho,li)
          #newgradf(initparm,y,M,ntimes,timeindex,udiff,zi,rho,li)
        })
        # proc.time() - time
      }else{
        
        cl <- parallel::makeCluster(ncores)
        parallel::clusterExport(cl,c("M","ylist","timelist","yceil","l","parm",
                                     "z","newrho","method",
                                     "newf","newgradf","grad_zipnegloglik_nocov_cont",
                                     "zipnegloglik_nocov_cont","J","libpath",
                                     "totalgroup","clusterindex","commonindex"),
                                envir=environment())
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
       
          fullindex <-  commonindex  #for z and l
          zi <- z[fullindex]
          li <- l[i,fullindex]
          initparm <- parm[i,]
          
          optim(par=initparm,fn=newf,gr=newgradf,
                M=M,y=y,ntimes=ntimes,timeindex=timeindex,udiff=udiff, 
                zi=zi,rho=newrho,li=li,
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
      
      #update z
      tempcommon <- J[jcommonindex,commonindex2]%*%t(parm[,commonindex2])
      z[commonindex] <- rowMeans(tempcommon) + colMeans(l[,commonindex])/rho
      newdiff <- sum((tempcommon-z[commonindex])^2)
       
      relchange <- newdiff / (1+olddiff)
      
      resid <- abs(newdiff-olddiff) / (1+olddiff)  
      resid_change <- c(resid_change, resid)
      #
      newnorm <- sum(z^2)
      primal_diff <- abs(sqrt(newnorm) - sqrt(oldnorm)) / (1+sqrt(oldnorm))
      primal_change <- c(primal_change, primal_diff)
      
      #update l
      for(i in 1:nsubj){
       
        fullindex <-  commonindex #for z and l
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
    return(list(working_parm=parm,
    #common_parm=z,
                threshold = c(primal_change[-1],dual_change[-1],nllk_change[-1]),
                nllk=nllk))
    
    ##################
    
  }else{ #some clustering some common
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
    zipnegloglik_nocov_cont <- ziphsmm::zipnegloglik_nocov_cont
    newf <- function(initparm,y,M,ntimes,timeindex,udiff,
                     zi,rho,li){
      part1 <- zipnegloglik_nocov_cont(initparm,M, y,ntimes,timeindex,udiff)
      diff <- J%*%initparm - zi
      part2 <- t(li)%*%diff
      part3 <- 0.5*rho*t(diff)%*%diff
      return(part1+part2+part3)
    }
    
    grad_zipnegloglik_nocov_cont <- ziphsmm::grad_zipnegloglik_nocov_cont
    newgradf <- function(initparm,y,M,ntimes,timeindex,udiff,
                         zi,rho,li){
      part1 <- grad_zipnegloglik_nocov_cont(initparm,M, y,ntimes,timeindex,udiff)
      part2 <- t(J)%*%li
      part3 <- rho * (t(J) %*% ( J%*%initparm - zi))
      return(part1+part2+part3)
    }
    
    #start iterations
    iteration <- 1
    nllk <- 0
    resid_change <- NULL
    dual_change <- NULL
    nllk_change <- NULL
    primal_change <- NULL
    
    
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
            
            #get subject-specific z and l
            for(kk in 1:totalgroup)
              if(i%in%group[[kk]]){gi <- kk}else{next}
            fullindex <-  c(clusterindex[[gi]],commonindex)  #for z and l
            zi <- z[fullindex]
            li <- l[i,fullindex]
            initparm <- parm[i,]
            
            optim(par=initparm,fn=newf,gr=newgradf,
                  M=M,y=y,ntimes=ntimes,timeindex=timeindex,udiff=udiff, 
                  zi=zi,rho=newrho,li=li,
                  method=method,...)
          #newf(initparm,y,M,ntimes,timeindex,udiff,zi,rho,li)
          #newgradf(initparm,y,M,ntimes,timeindex,udiff,zi,rho,li)
        })
        # proc.time() - time
      }else{
        
        cl <- parallel::makeCluster(ncores)
        parallel::clusterExport(cl,c("M","ylist","timelist","yceil","l","parm",
                                     "z","newrho","method","group",
                                     "newf","newgradf","grad_zipnegloglik_nocov_cont",
                                     "zipnegloglik_nocov_cont","J","libpath",
                                     "totalgroup","clusterindex","commonindex"),
                                     envir=environment())
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
            for(kk in 1:totalgroup)
              if(i%in%group[[kk]]){gi <- kk}else{next}
            fullindex <-  c(clusterindex[[gi]],commonindex)  #for z and l
            zi <- z[fullindex]
            li <- l[i,fullindex]
            initparm <- parm[i,]
          
          optim(par=initparm,fn=newf,gr=newgradf,
                M=M,y=y,ntimes=ntimes,timeindex=timeindex,udiff=udiff, 
                zi=zi,rho=newrho,li=li,
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
      newnorm <- sum(z^2)
      primal_diff <- abs(sqrt(newnorm) - sqrt(oldnorm)) / (1+sqrt(oldnorm))
      primal_change <- c(primal_change, primal_diff)
      
      relchange <- newdiff / (1+olddiff)
      resid <- abs(newdiff-olddiff) / (1+olddiff)  
      resid_change <- c(resid_change, resid)
      
      
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
      
      kkt_cur <- max(primal_diff, reldualchange, new_nllk_change) #newzdiff
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
    return(list(working_parm=parm,
    #common_parm=z,
                threshold=pmax(primal_change[-1],dual_change[-1],nllk_change[-1]),
                nllk=nllk))
  }
  
}



