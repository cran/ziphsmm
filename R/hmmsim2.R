
#######################################################
#' Simulate a hidden Markov series and its underlying states with covariates 
#' @param workparm working parameters
#' @param M number of latent states
#' @param n length of the simulated series
#' @param zeroindex a vector specifying whether a certain state is zero-inflated
#' @param prior_x matrix of covariates for generalized logit of prior probabilites (excluding the 
#' 1st probability). Default to NULL.
#' @param tpm_x matrix of covariates for transition probability matrix (excluding the 1st column).
#' Default to NULL.
#' @param emit_x matrix of covariates for the log poisson means. Default to NULL.
#' @param zeroinfl_x matrix of covariates for the nonzero structural zero proportions. Default to NULL.
#' @return a matrix with 1st column of simulated series and 
#' 2nd column of corresponding states
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' \dontrun{
#' priorparm <- 0
#' tpmparm <- c(0,-0.5,0.5,0,-0.2,0.8)
#' zeroindex <- c(1,0)
#' zeroparm <- c(0,-1,1)
#' emitparm <- c(2,0.5,-0.5,3,0.3,-0.2)
#' workparm <- c(priorparm,tpmparm,zeroparm,emitparm)
#' 
#' designx <- matrix(rnorm(2000),nrow=1000,ncol=2)
#' result <- hmmsim2(workparm,2,1000,zeroindex,tpm_x=designx,
#'       emit_x=designx,zeroinfl_x=designx)
#' 
#' y <- result$series
#' 
#' prior_init <- c(0.5,0.5)
#' emit_init <- c(10,30)
#' zero_init <- c(0.6,0)
#' omega <- matrix(c(0.9,0.1,0.2,0.8),2,2,byrow=TRUE)
#' 
#' 
#' fit <-  hmmfit(y,NULL,2,prior_init,omega,
#'      emit_init,zero_init, emit_x=designx,zeroinfl_x=designx,
#'      tpm_x=designx,hessian=FALSE,
#'      method="Nelder-Mead", control=list(maxit=2000,trace=1))
#' 
#' decode <- hmmviterbi2(y,NULL,2,fit$working_parameters,zero_init=c(1,0),
#'             emit_x=designx,zeroinfl_x=designx, tpm_x=designx,
#'             plot=TRUE, xlab="time", ylab="count",
#'             xlim=c(0,360),ylim=c(0,200))
#' 
#' sum(decode!=result$state)
#' }
#' 
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to fit a homogeneous hidden semi-markov model
hmmsim2 <- function(workparm, M, n, zeroindex, prior_x=NULL, tpm_x=NULL, emit_x=NULL, 
                   zeroinfl_x=NULL){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")

  #check if covariates
  tempmat <- matrix(0,nrow=n,ncol=1) #for NULL covariates
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
    stop("Need at least 1 covariates!")
  }else{
    #if there are covariates
    result <- hmm_cov_gen(workparm, M, n, ncovprior, covprior, ncovtpm, covtpm,
                      ncovzeroinfl, covzeroinfl, ncovemit, covemit, zeroindex)
    series <- result[,1]
    state <- result[,2]
    
    return(list(series=series,state=state))
    
  }
  
}
