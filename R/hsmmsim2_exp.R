
#######################################################
#' Simulate a hidden semi-Markov series and its underlying states with covariates
#' @param prior a vector of prior probabilities
#' @param dtrate a vector for the scale parameters in the base exponential
#' density for the latent state durations.
#' @param dtparm a matrix of coefficients for the accelerated failure time
#' model in each latent state
#' @param zeroparm a vector of regression coefficients for the structural
#' zero proportion in state 1
#' @param emitparm a matrix of regression coefficients for the Poisson
#' regression in each state
#' @param tpmparm a vector of coefficients for the multinomial logistic
#' regression in the transition probabilities
#' @param trunc a vector
#' @param M number of latent states
#' @param n length of the simulated series
#' @param dt_x if dt_dist is "nonparametric", then dt_x is the matrix of nonparametric 
#' state durataion probabilities. Otherwise, dt_x is matrix of covariates for the dwell time distribution 
#' parameters in log-series or shifted-poisson distributions.Default to NULL. 
#' @param tpm_x matrix of covariates for transition probability matrix (excluding the 1st column).
#' Default to NULL.
#' @param emit_x matrix of covariates for the log poisson means. Default to NULL.
#' @param zeroinfl_x matrix of covariates for the nonzero structural zero proportions. Default to NULL.
#' @return simulated series and corresponding states
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' 
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to fit a homogeneous hidden semi-markov model
hsmmsim2_exp <- function(prior,dtrate,dtparm,zeroparm,emitparm,tpmparm,
                         trunc, M, n, dt_x=NULL,
                         tpm_x=NULL, emit_x=NULL, zeroinfl_x=NULL){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  
  #check if covariates

  if(is.null(dt_x)) {
    covdt <- matrix(0,nrow=n,ncol=1)
  }else{
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
    covemit <- matrix(1,nrow=n,ncol=1)
  }else{
    covemit <- cbind(1,emit_x)
  }
  if(is.null(zeroinfl_x)){
    covzeroinfl <- matrix(1,nrow=n,ncol=1)
  }else{
    covzeroinfl <- cbind(1,zeroinfl_x)
  }
  
 #if there are covariates
  result <- hsmm_cov_gen_aft(prior, dtrate, dtparm, tpmparm,
                         zeroparm, emitparm, trunc,
                         M, n, ncovtpm, 
                        covdt, covtpm,covzeroinfl, covemit)
    
    series <- result[,1]
    state <- result[,2]
    
    return(list(series=series,state=state))
    
}



