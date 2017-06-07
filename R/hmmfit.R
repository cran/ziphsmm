
#######################################################
#' Estimate the parameters of a general zero-inflated Poisson hidden Markov model
#' by directly minimizing of the negative log-likelihood function using the
#' gradient descent algorithm.
#'
#' @param y observed time series values
#' @param ntimes a vector specifying the lengths of individual, 
#' i.e. independent, time series. If not specified, the responses are assumed to 
#' form a single time series, i.e. ntimes=length(y)
#' @param M number of latent states
#' @param prior_init a vector of initial values for prior probability for each state
#' @param tpm_init a matrix of initial values for transition probability matrix
#' @param emit_init a vector of initial values for the means for each poisson distribution
#' @param zero_init a vector of initial values for the structural zero proportions in each state
#' @param prior_x matrix of covariates for generalized logit of prior probabilites (excluding the 
#' 1st probability). Default to NULL.
#' @param tpm_x matrix of covariates for transition probability matrix (excluding the 1st column).
#' Default to NULL.
#' @param emit_x matrix of covariates for the log poisson means. Default to NULL.
#' @param zeroinfl_x matrix of covariates for the nonzero structural zero proportions. Default to NULL.
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
#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(10,40,70)
#' zero_init <- c(0.5,0,0)
#' omega <- matrix(c(0.5,0.3,0.2,0.4,0.3,0.3,0.2,0.4,0.4),3,3,byrow=TRUE)
#' result <- hmmsim(n=1000,M=3,prior=prior_init, tpm_parm=omega,
#'          emit_parm=emit_init,zeroprop=zero_init)
#' y <- result$series
#' fit <- hmmfit(y=y,M=3,prior_init=prior_init,tpm_init=omega,
#'      emit_init=emit_init,zero_init=zero_init,
#'      method="Nelder-Mead",hessian=TRUE,control=list(maxit=500,trace=1))
#' str(fit)
#' 
#' #variances for the 12 working parameters, which are the logit of prior probabilities
#' #for state 2 and state 3, the generalized logit of the transition probability matrix
#' #(tpm[1,2],tpm[1,3],tpm[2,2],tpm[2,3],tpm[3,2],tpm[3,3]), the logit of structural zero
#' #proportions for state 1, and log of poisson means for state 1, 2, and 3
#' #logit of non-zero zero proportions, and the log of poisson means
#' variance <- diag(solve(fit$obsinfo))
#' 
#' #with covariates
#' data(CAT)
#' y <- CAT$activity
#' x <- data.matrix(CAT$night)
#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(10,50,100)
#' zero_init <- c(0.5,0,0)
#' omega <- matrix(c(0.5,0.3,0.2,0.4,0.3,0.3,0.2,0.4,0.4),3,3,byrow=TRUE)
#' fit2 <-  hmmfit(y,rep(1440,3),3,prior_init,omega,
#'      emit_init,zero_init, emit_x=x,zeroinfl_x=x,hessian=FALSE,
#'      method="Nelder-Mead", control=list(maxit=500,trace=1))
#' fit2
#'      
#' \dontrun{
#' #two zero-inflated poissons
#' prior_init <- c(0.5,0.5)
#' emit_init <- c(10,50)
#' zero_init <- c(0.6,0.3)
#' omega <- matrix(c(0.9,0.1,0.2,0.8),2,2,byrow=TRUE)
#' result <- hmmsim(n=1000,M=2,prior=prior_init, tpm_parm=omega,
#'          emit_parm=emit_init,zeroprop=zero_init)
#' y <- result$series
#' fit <- hmmfit(y=y,M=2,prior_init=prior_init,tpm_init=omega,
#'      emit_init=emit_init,zero_init=zero_init,
#'      method="Nelder-Mead",hessian=FALSE,control=list(maxit=500,trace=1))
#' str(fit)
#' 
#' #four regular poissons
#' prior_init <- c(0.4,0.2,0.2,0.2)
#' emit_init <- c(10,40,70,100)
#' zero_init <- c(0,0,0,0)
#' omega <- matrix(c(0.3,0.3,0.2,0.2,0.4,0.2,0.3,0.1,
#'                   0.2,0.2,0.3,0.3,0.1,0.2,0.3,0.4),4,4,byrow=TRUE)
#' result <- hmmsim(n=1000,M=4,prior=prior_init, tpm_parm=omega,
#'          emit_parm=emit_init,zeroprop=zero_init)
#' y <- result$series
#' fit <- hmmfit(y=y,M=4,prior_init=prior_init,tpm_init=omega,
#'      emit_init=emit_init,zero_init=zero_init,
#'      method="Nelder-Mead",hessian=FALSE,control=list(maxit=500,trace=1))
#' str(fit)
#'}
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to fit a homogeneous hidden semi-markov model
hmmfit <- function(y, ntimes=NULL, M, prior_init, tpm_init,  emit_init, 
                  zero_init, prior_x=NULL, tpm_x=NULL, emit_x=NULL, 
                  zeroinfl_x=NULL,method="Nelder-Mead", hessian=FALSE, ...){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(prior_init)!=M | length(emit_init)!=M | length(zero_init)!=M |
     nrow(tpm_init)!= M | ncol(tpm_init)!=M) stop("The dimension of the initial value does not equal M!") 
  
  if(is.null(ntimes)) ntimes <- length(y)
  #check if covariates
  tempmat <- matrix(0,nrow=length(y),ncol=1) #for NULL covariates
  if(is.null(prior_x)) {
    ncovprior <- 0
    covprior <- tempmat
  }else{
    ncovprior <- ncol(prior_x)
    covprior <- prior_x
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
  
  #if no covariates
  if(ncovprior+ncovtpm+ncovemit+ncovzeroinfl==0){
    
    #obtain initial values for working parameters
    allparm <- rep(NA, M-1+M*(M-1)+sum(zero_init!=0)+M)
    lastindex <- 0
    allparm[(lastindex+1):(lastindex+M-1)] <- glogit(prior_init) #except 1st prior prob
    lastindex <- M-1
    for(j in 1:M){
      allparm[(lastindex+1):(lastindex+M-1)] <- glogit(tpm_init[j,])
      lastindex <- lastindex + M - 1
    }
    lastindex <- M*M-1
    for(i in 1:M) {
      if(zero_init[i]!=0) {
        allparm[lastindex+1] <- glogit(zero_init[i])
        lastindex <- lastindex + 1
      }
    }
    
    for(i in 1:M) allparm[lastindex+i] <- log(emit_init[i])
    
   
    #Maximum likelihood 
    parm <- optim(allparm,hmm_common_nocov_negloglik,M=M,ally=y,
                  ntimes=ntimes,zeroindex=zero_init,method=method,hessian=hessian,...)
    
    nllk <- parm$value
    aic <- nllk * 2 + 2 * (M*M+M)
    bic <- nllk * 2 + log(sum(ntimes)) * (M*M+M)
    workparm <- parm$par
    #use multivariate-delta method to get se's together
    obsinfo <- parm$hessian
    
    #retrieve the natural parameters from the working parameters
    lastindex <- 0
    prior <- ginvlogit(workparm[(lastindex+1):(lastindex+M-1)])
    lastindex <- M-1
    tpm <- matrix(0,M,M)
    for(j in 1:M){
      tpm[j,] <- ginvlogit(workparm[(lastindex+1):(lastindex+M-1)])
      lastindex <- lastindex+M-1
    }
    lastindex <- M*M-1
    
    zeroprop <- rep(NA, M)
    for(i in 1:M){
      if(zero_init[i]==0) zeroprop[i] <- 0
      else{
       zeroprop[i] <- exp(workparm[lastindex+1])/(1+exp(workparm[lastindex+1]))
       lastindex <- lastindex + 1
      }
    }
  
    emit_parm <- exp(workparm[(lastindex+1):(lastindex+M)])
 
  return(list(nllk=nllk,aic=aic,bic=bic,obsinfo=obsinfo,
              prior=prior,zeroprop=zeroprop,emit_parm=emit_parm,tpm=tpm))
  }else{
    #if there are covariates
    #obtain initial values for working parameters
    allparm <- rep(0, (M-1)*(ncovprior+1)+M*(M-1)*(ncovtpm+1)+
                     sum(zero_init!=0)*(ncovzeroinfl+1)+M*(ncovemit+1))
    lastindex <- 0
    allparm[seq(lastindex+1, lastindex+1+(M-2)*(ncovprior+1),length=M-1)] <- glogit(prior_init) #except 1st prior prob
    lastindex <- lastindex + (M-1)*(ncovprior+1)
    
    for(j in 1:M){
      allparm[seq(lastindex+1,lastindex+1+(M-2)*(ncovtpm+1),length=M-1)] <- 
        glogit(tpm_init[j,])
      lastindex <- lastindex + (M - 1)*(ncovtpm+1)
    }
    
    for(j in 1:M){
      if(zero_init[j]!=0){
        allparm[lastindex+1+(j-1)*(ncovzeroinfl+1)] <- glogit(zero_init[j])
        lastindex <- lastindex + ncovzeroinfl + 1
      }
    }
    
    for(i in 1:M) allparm[lastindex+1+(i-1)*(ncovemit+1)] <- log(emit_init[i])
    
    #hmm_cov_negloglik(allparm,3,y,"log","poisson",TRUE,
    #       1,x,0,x,0,x,1,x,1,x)
    #hmm_common_negloglik(allparm,3,y,trunc,c(1000,1000),"log","poisson",TRUE,
    #       1,x,0,x,0,x,1,x,1,x)
    #Maximum likelihood 
    parm <- optim(allparm,hmm_common_negloglik,M=M,ally=y,
                  ntimes=ntimes,zeroindex=zero_init,
                  ncolcovpi=ncovprior,allcovpi=covprior,ncolcovtrans=ncovtpm,
                  allcovtrans=covtpm,ncolcovp1=ncovzeroinfl,allcovp1=covzeroinfl,
                  ncolcovpois=ncovemit,allcovpois=covemit,
                  method=method,hessian=hessian,...) 

    #retrieve parameters
    nllk <- parm$value
    aic <- 2*nllk + 2*length(allparm)
    bic <- 2*nllk + log(sum(ntimes)) * length(allparm)
    result <- vector(mode="list", length=7)
    
    names(result) <- c("NLLK_AIC_BIC","observed_information",
                       "glogit_prior(except the 1st state)","glogit_tpm_parm(except the 1st column)",
                       "logit_zero_proportion(for nonzero ones)","log_emit_parm",
                       "working_parameters")
    result[[1]] <- c(nllk,aic,bic)
    #use multivariate-delta method to get se's together
    temp <- parm$hessian
    if(is.null(temp)) result[[2]]="None"
    else result[[2]] <- temp
    lastindex <- 0
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
    result[[3]] <- temp 
    
    lastindex <- lastindex + (M-1)*(ncovprior+1)
    temp <- matrix(0, nrow=M*(M-1), ncol=ncovtpm+1)
    rownames <- NULL
    for(i in 1:(M*(M-1))){
      rownames <- c(rownames,paste('glogit tpm_parm_',i,sep=""))
      temp[i,] <- parm$par[(lastindex+1+(i-1)*(ncovtpm+1)):(lastindex+i*(ncovtpm+1))]
    }
    colnames <- NULL
    for(j in 1:(ncovtpm+1)) colnames <- c(colnames,(paste('beta_',j-1,sep="")))
    rownames(temp) <- rownames
    colnames(temp) <- colnames
    result[[4]] <- temp
    
    lastindex <- lastindex + M*(M-1)*(ncovtpm+1)
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
    
    result[[7]] <- parm$par
    return(result)
  }
  
}
