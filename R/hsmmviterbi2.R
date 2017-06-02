
########################################################
#' Viterbi algorithm to decode the latent states in hidden semi-Markov models
#' with covariates
#' @param y the observed series to be decoded
#' @param ntimes vector specifying the lengths of individual, 
#' i.e. independent, time series. If not specified, the responses are assumed to 
#' form a single time series, i.e. ntimes=length(y)
#' @param M number of latent states
#' @param trunc a vector specifying the truncation at the maximum number of dwelling time 
#' in each state.
#' @param workparm working parameters. The first part is the logit of parameter p for 
#' each log-series dwell time distribution or the log of parameter theta for each 
#' shifted-poisson dwell time distribution. If dt_dist is "nonparametric", then the first
#' part is empty. The next part is the generalized logit of prior probabilities (except for 
#' the 1st state),the logit of each nonzero structural zero proportions, the log of each 
#' state-dependent poisson means, and the generalized logit of the transition probability 
#' matrix(except 1st column and the diagonal elements).
#' @param dt_dist dwell time distribution, can only be "log", "shiftedpoisson" or
#' "nonparametric". Default to "nonparametric".
#' @param zero_init a vector containing structural zero proportions in each state, e.g. set
#' zero_init[i] to be 0 if the i-th state is a regular poisson, and otherwise 1.
#' @param prior_x matrix of covariates for generalized logit of prior probabilites (excluding the 
#' 1st probability). Default to NULL.
#' @param tpm_x matrix of covariates for transition probability matrix (excluding the 1st column).
#' Default to NULL.
#' @param emit_x matrix of covariates for the log poisson means. Default to NULL.
#' @param zeroinfl_x matrix of covariates for the nonzero structural zero proportions. Default to NULL.
#' @param dt_x if dt_dist is "nonparametric", then dt_x is the matrix of nonparametric 
#' state durataion probabilities. Otherwise, dt_x is matrix of covariates for the dwell time distribution 
#' parameters in log-series or shifted-poisson distributions.Default to NULL.
#' @param plot whether a plot should be returned
#' @param xlim vector specifying the minimum and maximum on the x-axis in the plot. 
#' Default to NULL.
#' @param ylim vector specifying the minimum and maximum on the y-axis in the plot. 
#' Default to NULL.
#' @param ... further arguments to be passed to the plot() function 
#' @return decoded series of latent states
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' 
#' \dontrun{
#' data(CAT)
#' y <- CAT$activity
#' x <- data.matrix(CAT$night)
#' prior_init <- c(0.5,0.3,0.2)
#' dt_init <- c(0.9,0.6,0.3)
#' emit_init <- c(10,50,100)
#' zero_init <- c(0.5,0,0) #assuming only the 1st state has structural zero's
#' tpm_init <- matrix(c(0,0.3,0.7,0.4,0,0.6,0.5,0.5,0),3,3,byrow=TRUE)
#' trunc <- c(10,7,4)
#' fit2 <-  hsmmfit(y,rep(1440,3),3,trunc,prior_init,"log",dt_init,tpm_init,
#'      emit_init,zero_init,emit_x=x,zeroinfl_x=x,hessian=FALSE,
#'      method="Nelder-Mead", control=list(maxit=500,trace=1))
#' decode <- hsmmviterbi2(y,rep(1440,3),3,trunc,fit2$working_parameters,
#'             dt_dist="log", zero_init=c(1,0,0),
#'             emit_x=x,zeroinfl_x=x, plot=TRUE, xlab="time", ylab="count",
#'             xlim=c(0,360),ylim=c(0,200))
#' }
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hsmmviterbi2 <- function(y, ntimes=NULL, M, trunc, workparm, dt_dist="nonparametric",
                        zero_init,
                        prior_x=NULL, tpm_x=NULL, emit_x=NULL,zeroinfl_x=NULL,
                        dt_x=NULL,plot=FALSE, xlim=NULL, ylim=NULL, ...){
  if(!dt_dist%in%c("log","shiftpoisson","nonparametric"))
    stop("dt_dist can only be 'log' or 'shiftpoisson'!")
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(zero_init)!=M | length(trunc)!=M) 
    stop("The dimension of the zero_init or trunc does not equal M!") 
  
  if(is.null(ntimes)) ntimes <- length(y)
  state <- rep(NA, sum(ntimes))
  
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
  
  
  lastindex <- 0
  for(i in 1:length(ntimes)){
    if(dt_dist=="nonparametric"){
      state[(lastindex+1):(lastindex+ntimes[i])] <- hsmm_cov_viterbi_np(workparm,dt_x,
                                                                        M,trunc,y[(lastindex+1):(lastindex+ntimes[i])], 
                                                                        zero_init,
                                                                        ncovprior,covprior[(lastindex+1):(lastindex+ntimes[i]),,drop=FALSE],
                                                                        ncovtpm,covtpm[(lastindex+1):(lastindex+ntimes[i]),,drop=FALSE],
                                                                        ncovzeroinfl,covzeroinfl[(lastindex+1):(lastindex+ntimes[i]),,drop=FALSE],
                                                                        ncovemit,covemit[(lastindex+1):(lastindex+ntimes[i]),,drop=FALSE])
    }else{
      state[(lastindex+1):(lastindex+ntimes[i])] <- hsmm_cov_viterbi(workparm, M,
                                                                     y[(lastindex+1):(lastindex+ntimes[i])], trunc,zero_init,
                                                                     ncovdt,covdt[(lastindex+1):(lastindex+ntimes[i]),,drop=FALSE],
                                                                     ncovprior,covprior[(lastindex+1):(lastindex+ntimes[i]),,drop=FALSE],
                                                                     ncovtpm,covtpm[(lastindex+1):(lastindex+ntimes[i]),,drop=FALSE],
                                                                     ncovzeroinfl,covzeroinfl[(lastindex+1):(lastindex+ntimes[i]),,drop=FALSE],
                                                                     ncovemit,covemit[(lastindex+1):(lastindex+ntimes[i]),,drop=FALSE],
                                                                     dt_dist)
    }
    
    lastindex <- lastindex+ntimes[i]
  }
  
  if(plot==TRUE){
    xlimit <- rep(NA,2)
    ylimit <- rep(NA,2)
    if(is.null(xlim)){
      xlimit[1] <- 0
      xlimit[2] <- sum(ntimes)
    }else{xlimit <- xlim}
    if(is.null(ylim)){
      ylimit[1] <- 0
      ylimit[2] <- max(y) * 1.3
    }else{ylimit <- ylim}
    
    temp <- y[1:min(length(y),floor(xlimit[2]))]
    plot(temp, xlim=xlimit, ylim=ylimit, type="l",...)
    points(1:length(temp),rep(ylimit[1],length(temp)),
           pch=16,cex=0.8,col=state+1)
    legend <- NULL
    for(i in 1:M) legend <- c(legend,paste("state ",i,sep=""))
    legend("top",legend,pch=16,col=2:(M+1),bty="n",cex=0.9,
           ncol=M,xpd=TRUE)
  }    
  
  return(state)
}
