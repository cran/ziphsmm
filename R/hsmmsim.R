


######################################################################
#' Simulate a hidden semi-Markov series and its corresponding states
#' according to the specified parameters
#'
#' @param n length of the simulated series
#' @param M number of hidden states
#' @param prior a vector of prior probability for each state
#' @param dt_dist dwell time distribution, which should be "log" or
#' "shiftpoisson" or "nonparametric". Default to "nonparametric".
#' @param dt_parm a vector of dwell time distribution parameters for each state. If 
#' dt_dist is "log", then dt_parm is vector of p's; if dt_dist is "shiftpoisson", then
#' dt_parm is vector of theta's; if dt_dist is "nonparametric", then dt_parm is
#' a matrix whose i,j th element is the probability of staying in state i for duration j.
#' @param tpm_parm transition probability matrix, whose diagonal should be zero's.
#' @param emit_parm a vector containing means for each poisson distribution
#' @param zeroprop a vector containing structural zero proportions in each state 
#' @return simulated series and corresponding states
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' prior_init <- c(0.5,0.2,0.3)
#' dt_init <- c(0.8,0.5,0.2)
#' emit_init <- c(10,50,100)
#' zeroprop <- c(0.6,0.3,0.1)
#' omega <- matrix(c(0,0.3,0.7,0.4,0,0.6,0.5,0.5,0),3,3,byrow=TRUE)
#' sim1 <- hsmmsim(n=1000,M=3,prior=prior_init,dt_dist="log",
#'          dt_parm=dt_init, tpm_parm=omega,
#'          emit_parm=emit_init,zeroprop=zeroprop)
#' str(sim1)
#' 
#' prior_init <- c(0.5,0.5)
#' dt_init <- c(10,5)
#' emit_init <- c(10,30)
#' zeroprop <- c(0.5,0)
#' omega <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
#' sim2 <- hsmmsim(n=1000,M=2,prior=prior_init,dt_dist="shiftpoisson",
#'          dt_parm=dt_init, tpm_parm=omega,
#'          emit_parm=emit_init,zeroprop=zeroprop)
#' str(sim2)
#' hist(sim2$series,main="Histogram of observed values",xlab="observed values")
#' 
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
#main function to simulate hidden semi-markov model
hsmmsim <- function(n, M, prior, dt_dist="nonparametric", dt_parm, 
                     tpm_parm, emit_parm, zeroprop){
  #RcppExports.R to find the cpp function needed!!!
  if(!dt_dist%in%c("log","shiftpoisson","nonparametric")) 
    stop("dt_dist can only be 'log' or 'shiftpoisson' or 'nonparametric'!")
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(prior)!=M | length(emit_parm)!=M | length(zeroprop)!=M |
     nrow(tpm_parm)!= M | ncol(tpm_parm)!=M) stop("The dimension of the initial value does not equal M!") 
  
  if(dt_dist=="nonparametric"){
    result <- hsmm_gen_np(dim=n,M=M,pi=prior,theta=emit_parm,
                          zeroprop=zeroprop,omega=tpm_parm,dt=dt_parm)
  }
  else{
    result <- hsmm_gen(dim=n, M=M, pi=prior, theta=emit_parm,
                       zeroprop=zeroprop, omega=tpm_parm, p=dt_parm,
                       dt_dist=dt_dist)
  }
    series <- result[,1]
    state <- result[,2]
  
  return(list(series=series,state=state))
}
