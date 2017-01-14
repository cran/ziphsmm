
######################################################################
#' Simulate a hidden Markov series and its underlying states
#' with zero-inflated emission distributions
#'
#' @param n length of the simulated series
#' @param M number of hidden states
#' @param prior a vector of prior probability for each state
#' @param tpm_parm transition probability matrix
#' @param emit_parm a vector containing means for each poisson distribution
#' @param zeroprop a vector containing structural zero proportions in each state
#' @return a matrix with 1st column of simulated series and 
#' 2nd column of corresponding states
#' #' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(10,50,100)
#' zeroprop <- c(0.5,0,0)
#' omega <- matrix(c(0.5,0.3,0.2,0.4,0.3,0.3,0.2,0.4,0.4),3,3,byrow=TRUE)
#' result <- hmmsim(n=1000,M=3,prior=prior_init, tpm_parm=omega,
#'          emit_parm=emit_init,zeroprop=zeroprop)
#' str(result)
#' 
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to simulate hidden semi-markov model
hmmsim <- function(n, M, prior, tpm_parm, emit_parm, zeroprop){
  #RcppExports.R to find the cpp function needed!!!
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(prior)!=M | length(emit_parm)!=M | length(zeroprop)!=M |
     nrow(tpm_parm)!= M | ncol(tpm_parm)!=M) stop("The dimension of the initial value does not equal M!") 
  
    result <- hmm_gen(dim=n, M=M, ntimes=1, pi=prior, theta=emit_parm,
                       zeroprop=zeroprop, a=tpm_parm)
    series <- result[,1]
    state <- result[,2]
  
  return(list(series=series,state=state))
}
