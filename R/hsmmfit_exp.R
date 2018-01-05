
#######################################################
#' Simulate a hidden semi-Markov series and its underlying states with covariates 
#' where the latent state distributions have accelerated failure time structure
#' whose base densities are exponential
#' @param y observed time series values
#' @param M number of latent states
#' @param trunc a vector specifying truncation at the maximum number of dwelling time in 
#' each state.
#' @param dtrate a vector for the scale parameters in the base exponential
#' density for the latent state durations.
#' @param dtparm a matrix of coefficients for the accelerated failure time
#' model in each latent state
#' @param prior a vector of prior probabilities
#' @param zeroparm a vector of regression coefficients for the structural
#' zero proportion in state 1
#' @param emitparm a matrix of regression coefficients for the Poisson
#' regression in each state
#' @param tpmparm a vector of coefficients for the multinomial logistic
#' regression in the transition probabilities
#' @param dt_x a matrix of covariates for the latent state durations
#' @param zeroinfl_x a matrix of covariates for the zero proportion
#' @param emit_x a matrix of covariates for the Poisson means 
#' @param tpm_x a matrix of covariates for the transition
#' @param yceil a scalar defining the ceiling of y, above which the values will be
#' truncated. Default to NULL.
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
#' 
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to fit a homogeneous hidden semi-markov model
hsmmfit_exp <- function(y, M, trunc, dtrate, dtparm, prior,
                        zeroparm, emitparm, tpmparm,
                        dt_x, zeroinfl_x, emit_x, tpm_x,
                        yceil=NULL, method="Nelder-Mead", hessian=FALSE, ...){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(!is.null(yceil)) y <- ifelse(y>yceil, yceil, y)
  n <- length(y)
  #check if covariates
  if(is.null(dt_x)) {
    ncovdt <- 0
    covdt <- matrix(0,nrow=n,ncol=1)
  }else{
    ncovdt <- ncol(dt_x)
    covdt <- dt_x
  }
  
  #tpm_x, emit_x, zeroinfl_x
  if(is.null(tpm_x)){
    covtpm <- matrix(0,nrow=n,ncol=1)
    ncovtpm <- 0
  }else{
    covtpm <- tpm_x
    ncovtpm <- ncol(tpm_x)
  }
  if(is.null(emit_x)){
    ncovemit <- 0
    covemit <- matrix(1,nrow=n,ncol=1)
  }else{
    ncovemit <- ncol(emit_x)
    covemit <- cbind(1,emit_x)
  }
  if(is.null(zeroinfl_x)){
    ncovzeroinfl <- 0
    covzeroinfl <- matrix(1,nrow=n,ncol=1)
  }else{
    ncovzeroinfl <- ncol(zeroinfl_x)
    covzeroinfl <- cbind(1,zeroinfl_x)
  }
  
  if(ncovtpm==0 & ncovdt==0) {
    print("No covariates for state durations and transitions!")
    print("Use fasthsmmfit() to achieve better efficiency!")
  }else{ 
    
    if(M>2 & ncovdt>0) 
      allparm <- c(log(dtrate),dtparm, glogit(prior),
                    zeroparm, c(t(emitparm)),tpmparm)
    if(M>2 & ncovdt==0)
      allparm <- c(log(dtrate),        glogit(prior),
                   zeroparm, c(t(emitparm)),tpmparm)
    if(M==2 & ncovdt>0)
      allparm <- c(log(dtrate),dtparm, glogit(prior),
                   zeroparm, c(t(emitparm)))
    
    tempresult <- optim(allparm, hsmm_cov_nllk_aft, y=y, trunc=trunc,
                        M=M, ncolcovomega=ncovtpm, ncolcovp=ncovdt,
                        ncolcovp1=ncovzeroinfl, ncolcovpois=ncovemit,
                        covp=covdt,covomega=covtpm,covp1=covzeroinfl,
                        covpois=covemit,method=method,
                        hessian=hessian,...)
    fitparm <- tempresult$par
    negloglik <- tempresult$value
    natparm <- retrieve_hsmm_aft(fitparm,trunc,M,ncovtpm,ncovdt,
                                 ncovzeroinfl,ncovemit)
    
    return(list(prior=natparm$delta,dtrate=natparm$lden,
                dtparm=natparm$dtparm,tpmparm=natparm$tpmparm,
                zeroparm=natparm$zeroparm,emitparm=natparm$emitparm,
                working_parm=fitparm,negloglik=negloglik))
  }
  
}



