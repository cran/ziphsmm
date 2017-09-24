######################################################################
#' Simulate a hidden Markov series and its underlying states
#' with zero-inflated emission distributions
#'
#' @param n length of the simulated series
#' @param M number of hidden states
#' @param prior a vector of prior probability for each state
#' @param tpm_parm transition rate matrix
#' @param emit_parm a vector containing means for each poisson distribution
#' @param zeroprop a vector containing structural zero proportions in each state
#' @param timeindex a vector containing the time points
#' @return simulated series and corresponding states
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(10,40,70)
#' zero_init <- c(0.5,0,0)
#' omega <- matrix(c(-0.3,0.2,0.1,0.1,-0.2,0.1,0.2,0.2,-0.4),3,3,byrow=TRUE)
#' timeindex <- rep(1,1000)
#' for(i in 2:1000) timeindex[i] <- timeindex[i-1] + sample(1:3,1)

#' result <- hmmsim.cont(n=1000,M=3,prior=prior_init, tpm_parm=omega,
#'           emit_parm=emit_init,zeroprop=zero_init,timeindex=timeindex)
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hmmsim.cont <- function(n, M, prior, tpm_parm, emit_parm, zeroprop, timeindex){
  #RcppExports.R to find the cpp function needed!!!
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(prior)!=M | length(emit_parm)!=M | length(zeroprop)!=M |
     nrow(tpm_parm)!= M | ncol(tpm_parm)!=M) stop("The dimension of the initial value does not equal M!") 
  
  result <- hmm_gen_cont(dim=n, M=M, pi=prior, theta=emit_parm,
                    zeroprop=zeroprop, a=tpm_parm, timeindex=timeindex)
  series <- result[,1]
  state <- result[,2]
  
  return(list(series=series,state=state))
}

