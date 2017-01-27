
#######################################################
#' Simulate a hidden semi-Markov series and its underlying states with covariates 
#' @param workparm working parameters,which are the logit of parameter p for 
#' each log-series dwell time distribution or the log of parameter theta for each 
#' shifted-poisson dwell time distribution, the generalized logit of prior probabilities (except for the
#' 1st state),the logit of each nonzero structural zero proportions, the log of each 
#' state-dependent poisson means, and the generalized logit of the transition probability 
#' matrix(except 1st column and the diagonal elements)
#' @param M number of latent states
#' @param n length of the simulated series
#' @param zeroindex a vector specifying whether a certain state is zero-inflated
#' @param dt_dist dwell time distribution, can only be "log", "geometric",
#' or "shiftedpoisson"
#' @param dt_x matrix of covariates for the dwell time distribution parameters 
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
#' #example 1
#' 
#' dtparm <- c(6,2)
#' priorparm <- 0
#' zeroindex <- c(1,0)
#' zeroparm <- c(0,-1,1)
#' emitparm <- c(2,0.5,-0.5,4,0.3,-0.2)
#' workparm <- c(dtparm,priorparm,zeroparm,emitparm)
#' 
#' designx <- matrix(rnorm(4000),nrow=2000,ncol=2)
#' result <- hsmmsim2(workparm,2,2000,zeroindex,"shiftpoisson",
#'       emit_x=designx,zeroinfl_x=designx)
#' 
#' y <- result$series
#' 
#' dt_init <- c(8,3)
#' prior_init <- c(0.5,0.5)
#' emit_init <- c(10,50)
#' trunc <- c(13,8)
#' zeroprop <- c(0.8,0)
#' omega <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
#' 
#' 
#' fit1 <- hsmmfit(y=y,M=2,trunc=trunc,prior_init=prior_init,dt_dist="shiftpoisson",
#'      dt_init=dt_init,emit_x=designx,zeroinfl_x=designx,
#'      tpm_init=omega,emit_init=emit_init,zero_init=zeroprop,
#'      method="Nelder-Mead",hessian=FALSE,control=list(maxit=2000,trace=1))
#' 
#' decode <- hsmmviterbi2(y,NULL,2,trunc,fit1$working_parameters,
#'             dt_dist="log", zero_init=zeroprop,
#'             emit_x=designx,zeroinfl_x=designx, plot=TRUE, xlab="time", ylab="count",
#'             xlim=c(0,360),ylim=c(0,200))
#' sum(decode!=result$state)
#' }
#' 
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to fit a homogeneous hidden semi-markov model
hsmmsim2 <- function(workparm, M, n, zeroindex, dt_dist,  dt_x=NULL,prior_x=NULL, 
                     tpm_x=NULL, emit_x=NULL, zeroinfl_x=NULL){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(!dt_dist%in%c("log","shiftpoisson")) stop("dt_dist can only be 'log' or 'shiftpoisson'!")
  
  #check if covariates
  tempmat <- matrix(0,nrow=n,ncol=1) #for NULL covariates
  if(is.null(dt_x)) {
    ncovdt <- 0
    covdt <- tempmat
  }else{
    ncovdt <- ncol(dt_x)
    covdt <- dt_x
  }
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
  if(ncovdt+ncovprior+ncovtpm+ncovemit+ncovzeroinfl==0){
    stop("Need at least 1 covariates!")
  }else{
    #if there are covariates
    result <- hsmm_cov_gen(workparm, M, n, dt_dist, zeroindex, 
                           ncovdt, covdt, ncovprior, covprior, ncovtpm, covtpm,
                          ncovzeroinfl, covzeroinfl, ncovemit, covemit)
    series <- result[,1]
    state <- result[,2]
    
    return(list(series=series,state=state))
    
  }
  
}
