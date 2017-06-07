
#######################################################
#' Estimate the parameters of a general zero-inflated Poisson hidden semi-Markov model
#' by directly minimizing of the negative log-likelihood function using the gradient
#' descent algorithm.
#'
#' @param y observed time series values
#' @param ntimes A vector specifying the lengths of individual, 
#' i.e. independent, time series. If not specified, the responses are assumed to 
#' form a single time series, i.e. ntimes=length(data)
#' @param M number of hidden states
#' @param trunc a vector specifying truncation at the maximum number of dwelling time in each state. 
#' The higher the truncation, the more accurate the approximation but also the more 
#' computationally expensive.
#' @param prior_init a vector of initial value for prior probability for each state
#' @param dt_dist dwell time distribution, can only be "log", "geometric",
#' or "shiftedpoisson"
#' @param dt_init a vector of initial value for the parameter in each dwell time distribution, which 
#' should be a vector of p's for dt_dist == "log" and a vector of theta's for 
#' dt_dist=="shiftpoisson"
#' @param tpm_init a matrix of initial values for the transition probability matrix, whose diagonal
#' elements should be zero's
#' @param emit_init a vector initial value for the vector containing means for each poisson distribution
#' @param zero_init a vector initial value for the vector containing structural zero proportions in each state
#' @param prior_x matrix of covariates for generalized logit of prior probabilites (excluding the 
#' 1st probability). Default to NULL.
#' @param dt_x matrix of covariates for the dwell time distribution parameters 
#' @param tpm_x matrix of covariates for transition probability matrix (excluding the 1st column).
#' Default to NULL.
#' @param emit_x matrix of covariates for the log poisson means. Default to NULL.
#' @param zeroinfl_x matrix of covariates for the nonzero structural zero proportions. Default to NULL.
#' @param method method to be used for direct numeric optimization. See details in
#' the help page for optim() function. Default to Nelder-Mead.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' Note that the hessian is for the working parameters, which are the logit of parameter p for 
#' each log-series dwell time distribution or the log of parameter theta for each 
#' shifted-poisson dwell time distribution, the generalized logit of prior probabilities (except for the
#' 1st state),the logit of each nonzero structural zero proportions, the log of each 
#' state-dependent poisson means, and the generalized logit of the transition probability 
#' matrix(except 1st column and the diagonal elements)
#' @param ... Further arguments passed on to the optimization methods
#' @return simulated series and corresponding states
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' 
#' #2 zero-inflated poissons
#' prior_init <- c(0.5,0.5)
#' emit_init <- c(10,30)
#' dt_init <- c(10,6)
#' trunc <- c(20,10)
#' zeroprop <- c(0.5,0.3)
#' omega <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
#' sim2 <- hsmmsim(n=1000,M=2,prior=prior_init,dt_dist="shiftpoisson",
#'          dt_parm=dt_init, tpm_parm=omega,
#'          emit_parm=emit_init,zeroprop=zeroprop)
#' str(sim2)
#' y <- sim2$series
#' fit2 <- hsmmfit(y=y,M=2,trunc=trunc,prior_init=prior_init,dt_dist="shiftpoisson",
#'      dt_init=dt_init,
#'      tpm_init=omega,emit_init=emit_init,zero_init=zeroprop,
#'      method="Nelder-Mead",hessian=FALSE,control=list(maxit=500,trace=1))
#' str(fit2)
#' 
#' 
#' \dontrun{
#' #1 zero-inflated poisson and 3 regular poissons
#' prior_init <- c(0.5,0.2,0.2,0.1)
#' dt_init <- c(0.8,0.7,0.6,0.5)
#' emit_init <- c(10,30,70,130)
#' trunc <- c(10,10,10,10)
#' zeroprop <- c(0.6,0,0,0)  #only the 1st-state is zero-inflated
#' omega <- matrix(c(0,0.5,0.3,0.2,0.4,0,0.4,0.2,
#'                   0.2,0.6,0,0.2,0.1,0.1,0.8,0),4,4,byrow=TRUE)
#' sim1 <- hsmmsim(n=2000,M=4,prior=prior_init,dt_dist="log",
#'          dt_parm=dt_init, tpm_parm=omega,
#'          emit_parm=emit_init,zeroprop=zeroprop)
#' str(sim1)
#' y <- sim1$series
#' fit <- hsmmfit(y=y,M=4,trunc=trunc,prior_init=prior_init,dt_dist="log",dt_init=dt_init,
#'      tpm_init=omega,emit_init=emit_init,zero_init=zeroprop,
#'      method="Nelder-Mead",hessian=TRUE,control=list(maxit=500,trace=1))
#' str(fit)
#' 
#' #variances for the 20 working parameters, which are the logit of parameter p for 
#' #the 4 log-series dwell time distributions, the generalized logit of prior probabilities 
#' #for state 2,3,4, the logit of each nonzero structural zero proportions in state 1, 
#' #the log of 4 state-dependent poisson means, and the generalized logit of the 
#' #transition probability matrix(which are tpm[1,3],tpm[1,4], tpm[2,3],tpm[2,4],
#' #tpm[3,2],tpm[3,4],tpm[4,2],tpm[4,3])
#' variance <- diag(solve(fit$obsinfo)) 
#' 
#' 
#' #1 zero-inflated poisson and 2 poissons with covariates
#' data(CAT)
#' y <- CAT$activity
#' x <- data.matrix(CAT$night)
#' prior_init <- c(0.5,0.3,0.2)
#' dt_init <- c(0.9,0.6,0.3)
#' emit_init <- c(10,20,30)
#' zero_init <- c(0.5,0,0) #assuming only the 1st state has structural zero's
#' tpm_init <- matrix(c(0,0.3,0.7,0.4,0,0.6,0.5,0.5,0),3,3,byrow=TRUE)
#' trunc <- c(10,7,4)
#' fit2 <-  hsmmfit(y,rep(1440,3),3,trunc,prior_init,"log",dt_init,tpm_init,
#'      emit_init,zero_init,emit_x=x,zeroinfl_x=x,hessian=FALSE,
#'      method="Nelder-Mead", control=list(maxit=500,trace=1))
#' fit2
#' 
#' #another example with covariates for 2 zero-inflated poissons
#' data(CAT)
#' y <- CAT$activity
#' x <- data.matrix(CAT$night)
#' prior_init <- c(0.5,0.5)
#' dt_init <- c(10,5)
#' emit_init <- c(10, 30)
#' zero_init <- c(0.5,0.2)
#' tpm_init <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
#' trunc <- c(10,5)
#' fit <-  hsmmfit(y,NULL,2,trunc,prior_init,"shiftpoisson",dt_init,tpm_init,
#'      emit_init,zero_init,dt_x=x,emit_x=x,zeroinfl_x=x,tpm_x=x,hessian=FALSE,
#'      method="Nelder-Mead", control=list(maxit=500,trace=1))
#' fit
#' }
#' 
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to fit a homogeneous hidden semi-markov model
hsmmfit <- function(y, ntimes=NULL, M, trunc, prior_init, dt_dist, dt_init, 
                   tpm_init, emit_init, zero_init,
                   prior_x=NULL, dt_x=NULL, tpm_x=NULL, emit_x=NULL, 
                   zeroinfl_x=NULL, method="Nelder-Mead", hessian=FALSE, ...){
  if(!dt_dist%in%c("log","shiftpoisson")) stop("dt_dist can only be 'log' or 'shiftpoisson'!")
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(prior_init)!=M | length(emit_init)!=M | length(zero_init)!=M |
     nrow(tpm_init)!= M | ncol(tpm_init)!=M | length(trunc)!=M |
     length(dt_init) != M) stop("The dimension of the initial value does not equal M!") 
  
  if(is.null(ntimes)) ntimes <- length(y)
  #check if there are covariates
  tempmat <- matrix(0,nrow=length(y),ncol=1) #for NULL covariates
  if(is.null(prior_x)) {
    ncovprior <- 0
    covprior <- tempmat
  }else{
    ncovprior <- ncol(prior_x)
    covprior <- prior_x
  }
  if(is.null(dt_x)) {
    ncovdt <- 0
    covdt <- tempmat
  }else{
    ncovdt <- ncol(dt_x)
    covdt <- dt_x
  }
  #tpm_x, emit_x, zeroinfl_x
  if(is.null(tpm_x)){
    ncovtpm <- 0
    covtpm <- tempmat
  }else{
    ncovtpm <- ncol(tpm_x)
    covtpm <- tpm_x
  }
  if(is.null(emit_x)){
    ncovemit <- 0
    covemit <- ncol(emit_x)
  }else{
    ncovemit <- ncol(emit_x)
    covemit <- emit_x
  }
  if(is.null(zeroinfl_x)){
    ncovzeroinfl <- 0
    covzeroinfl <- tempmat
  }else{
    ncovzeroinfl <- ncol(zeroinfl_x)
    covzeroinfl <- zeroinfl_x
  }
  if(M==2) ncovtpm=0
  
  #if no covariates
  if(ncovprior+ncovtpm+ncovemit+ncovzeroinfl==0){
    #obtain initial values for working parameters
    allparm <- rep(NA, M+M-1+sum(zero_init!=0)+M+M*(M-2))
    if(dt_dist=="log"){
       for(i in 1:M) allparm[i] <- glogit(dt_init[i])
    }else{
      for(i in 1:M) allparm[i] <- log(dt_init[i])
    }
    
    lastindex <- M
    allparm[(lastindex+1):(lastindex+M-1)] <- glogit(prior_init) #except 1st prior prob
    lastindex <- M+M-1
    
    for(i in 1:M) {
      if(zero_init[i]!=0) {
        allparm[lastindex+1] <- glogit(zero_init[i])
        lastindex <- lastindex + 1
      }
    }
    
    for(i in 1:M) allparm[lastindex+i] <- log(emit_init[i])
    lastindex <- lastindex + M
    
    if(M>2){
      for(j in 1:M){
        allparm[(lastindex+1):(lastindex+M-2)] <- glogit(tpm_init[j,-j])
        lastindex <- lastindex + M - 2
      }
    }
    
    #hsmm_nllk(allparm,3,trunc,y,"log","poisson",TRUE)
    #hsmm_common_nocov_nllk(allparm,3,y,trunc,c(500,500),"log","poisson",TRUE)
    #Maximum likelihood 
    parm <- optim(allparm,hsmm_common_nocov_nllk,M=M,ally=y,trunc=trunc,
                  ntimes=ntimes,dt_dist=dt_dist,zeroprop=zero_init,
                  method=method,hessian=hessian,...)
    
    nllk <- parm$value
    aic <- nllk * 2 + 2 * (M*M+M)
    bic <- nllk * 2 + log(sum(ntimes)) * (M*M+M)
    workparm <- parm$par
    #use multivariate-delta method to get se's together
    obsinfo <- parm$hessian
    
    #retrieve the natural parameters from the working parameters
    if(dt_dist=="log"){    
      dt_parm <- exp(workparm[1:M])/(1+exp(workparm[1:M]))
    }else{
      dt_parm <- exp(workparm[1:M])
    }
    
    lastindex <- M
    prior <- ginvlogit(workparm[(lastindex+1):(lastindex+M-1)])
    lastindex <- 2*M-1
    
    zeroprop <- rep(NA, M)
    for(i in 1:M){
      if(zero_init[i]==0) zeroprop[i] <- 0
      else{
        zeroprop[i] <- exp(workparm[lastindex+1])/(1+exp(workparm[lastindex+1]))
        lastindex <- lastindex + 1
      }
    }
    
    emit_parm <- exp(workparm[(lastindex+1):(lastindex+M)])
    lastindex <- lastindex + M
    
    if(M==2) {
      tpm <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
    }else{
      tpm <- matrix(0,M,M)
      for(j in 1:M){
        tpm[j,-j] <- ginvlogit(workparm[(lastindex+1):(lastindex+M-2)])
        lastindex <- lastindex+M-2
      }
    }
    return(list(nllk=nllk,aic=aic,bic=bic,obsinfo=obsinfo,dt_parm=dt_parm,
                prior=prior,zeroprop=zeroprop,emit_parm=emit_parm,tpm=tpm))
  }else{
    #if there are covariates
    if(M>2){
      allparm <- rep(0, M*(ncovdt+1)+(M-1)*(ncovprior+1)+sum(zero_init!=0)*(ncovzeroinfl+1)+
                       M*(ncovemit+1)+M*(M-2)*(ncovtpm+1))
    }else{
      allparm <- rep(0, M*(ncovdt+1)+(M-1)*(ncovprior+1)+sum(zero_init!=0)*(ncovzeroinfl+1)+
                       M*(ncovemit+1))
    }
    
    for(i in 1:M){
      if(dt_dist=="log"){
        allparm[(i-1)*(ncovdt+1)+1] <- glogit(dt_init[i])
      }else{
        allparm[(i-1)*(ncovdt+1)+1] <- log(dt_init[i])
      }
    }
    
    lastindex <- M*(ncovdt+1)
    
    allparm[seq(lastindex+1, lastindex+1+(M-2)*(ncovprior+1),length=M-1)] <- glogit(prior_init) #except 1st prior prob
    
    lastindex <- lastindex + (M-1)*(ncovprior+1)
    
    for(j in 1:M){
      if(zero_init[j]!=0){
        allparm[lastindex+1+(j-1)*(ncovzeroinfl+1)] <- glogit(zero_init[j])
        lastindex <- lastindex + ncovzeroinfl + 1
      }
    }
    
    for(i in 1:M) allparm[lastindex+1+(i-1)*(ncovemit+1)] <- log(emit_init[i])
    lastindex <- lastindex + M*(ncovemit+1)
    
    if(M>2){
      for(j in 1:M){
        allparm[seq(lastindex+1,lastindex+1+(M-3)*(ncovtpm+1),length=M-2)] <- 
          glogit(tpm_init[j,-j])
        lastindex <- lastindex + (M - 2)*(ncovtpm+1)
      }
    }
    
    #optimization
    parm <- optim(allparm,hsmm_common_negloglik,M=M,ally=y,trunc=trunc,
                  ntimes=ntimes,dt_dist=dt_dist,zeroindex=zero_init,
                  ncolcovp=ncovdt,allcovp=covdt,
                  ncolcovpi=ncovprior,allcovpi=covprior,ncolcovomega=ncovtpm,
                  allcovomega=covtpm,ncolcovp1=ncovzeroinfl,allcovp1=covzeroinfl,
                  ncolcovpois=ncovemit,allcovpois=covemit, 
                  method=method,hessian=hessian,...)
    
    #retrieve parameters
    nllk <- parm$value
    aic <- 2*nllk + 2*length(allparm)
    bic <- 2*nllk + log(sum(ntimes)) * length(allparm)
    result <- vector(mode="list", length=8)
    
    names(result) <- c("NLLK_AIC_BIC","observed_information",
                       "logit_dt_parm","glogit_prior_parm",
                       "logit_zero_proportion",
                       "log_emit_parm","glogit_tpm_parm",
                       "working_parameters")
    result[[1]] <- c(nllk,aic,bic)
    #use multivariate-delta method to get se's together
    temp <- parm$hessian
    if(is.null(temp)) result[[2]]="None"
    else result[[2]] <- temp
    lastindex <- 0
    
    temp <- matrix(0,nrow=M,ncol=ncovdt+1)
    rownames <- NULL
    for(i in 1:M) {
      rownames <- c(rownames,paste("logit dt_parm_",i,sep=""))
      temp[i,] <- parm$par[((i-1)*(ncovdt+1)+1):(i*(ncovdt+1))]
    }
    colnames <- NULL
    for(j in 1:(ncovdt+1)) colnames <- c(colnames,paste("beta_",j-1,sep=""))
    rownames(temp) <- rownames
    colnames(temp) <- colnames
    result[[3]] <- temp 
    
    lastindex <- M*(ncovdt+1)
    temp <- matrix(0,nrow=M-1,ncol=ncovprior+1)
    rownames <- NULL
    for(i in 1:(M-1)){
      rownames <- c(rownames,paste('glogit prior_parm_',i,sep=""))
      temp[i,] <- parm$par[(lastindex+1+(i-1)*(ncovprior+1)):(lastindex+i*(ncovprior+1))]
    }
    colnames <- NULL
    for(j in 1:(ncovprior+1)) colnames <- c(colnames,paste('beta_',j-1,sep=""))
    rownames(temp) <- rownames
    colnames(temp) <- colnames
    result[[4]] <- temp 
    
    lastindex <- lastindex + (M-1)*(ncovprior+1)
    temp <- matrix(0,nrow=sum(zero_init!=0),ncol=ncovzeroinfl+1)
    rownames <- NULL
    for(i in 1:(sum(zero_init!=0))){
      rownames <- c(rownames,paste("logit zero proportion_",i,sep=""))
      temp[i,] <- parm$par[(lastindex+1):(lastindex+ncovzeroinfl+1)]
      lastindex <- lastindex + ncovzeroinfl + 1
    }
    colnames <- NULL
    for(j in 1:(ncovzeroinfl+1)) colnames <- c(colnames,paste('beta_',j-1,sep=""))
    colnames(temp) <- colnames
    rownames(temp) <- rownames
    result[[5]] <- temp
    
    temp <- matrix(0,nrow=M, ncol=ncovemit+1)
    rownames <- NULL
    for(i in 1:M){
      rownames <- c(rownames,paste('log emit_parm_',i,sep=""))
      temp[i,] <- parm$par[(lastindex+1+(i-1)*(ncovemit+1)):(lastindex+i*(ncovemit+1))]
    }
    colnames <- NULL
    for(j in 1:(ncovemit+1)) colnames <- c(colnames,paste('beta_',j-1,sep=""))
    rownames(temp) <- rownames
    colnames(temp) <- colnames
    result[[6]] <- temp
    
    lastindex <- lastindex + M*(ncovemit + 1)
    if(M==2){
      result[[7]] <- "when M=2, the tpm is always a matrix with diagonal 1 and off-diagonal 0"
    }else{
      temp <- matrix(0, nrow=M*(M-2), ncol=ncovtpm+1)
      rownames <- NULL
      for(i in 1:(M*(M-2))){
        rownames <- c(rownames,paste('glogit tpm_parm_',i,sep=""))
        temp[i,] <- parm$par[(lastindex+1+(i-1)*(ncovtpm+1)):(lastindex+i*(ncovtpm+1))]
      }
      colnames <- NULL
      for(j in 1:(ncovtpm+1)) colnames <- c(colnames,(paste('beta_',j-1,sep="")))
      rownames(temp) <- rownames
      colnames(temp) <- colnames
      result[[7]] <- temp
    }
    
    result[[8]] <- parm$par
    return(result)
  }
  
  
}
