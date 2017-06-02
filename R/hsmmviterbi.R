
########################################################
#' Viterbi algorithm to decode the latent states for hidden semi-Markov models
#' @param y the observed series to be decoded
#' @param ntimes vector specifying the lengths of individual, 
#' i.e. independent, time series. If not specified, the responses are assumed to 
#' form a single time series, i.e. ntimes=length(y)
#' @param M number of latent states
#' @param trunc a vector specifying truncation at the maximum number of dwelling time in 
#' each state. 
#' @param prior a vector of prior probability values
#' @param dt_dist dwell time distribution, can only be "log", "shiftedpoisson", or
#' "nonparametric". Default to "nonparametric".
#' @param dt_parm a vector of parameter values in each dwell time distribution, which 
#' should be a vector of p's for dt_dist == "log", or a vector of theta's for 
#' dt_dist=="shiftpoisson", or a matrix whose a matrix whose i,j th element is the probability 
#' of staying in state i for duration j for dt_dist == "nonparametric".
#' @param tpm_parm transition probability matrix
#' @param emit_parm a vector containing means for each poisson distribution
#' @param zero_init a vector containing structural zero proportions in each state
#' @param plot whether a plot should be returned
#' @param xlim vector specifying the minimum and maximum on the x-axis in the plot. 
#' Default to NULL.
#' @param ylim vector specifying the minimum and maximum on the y-axis in the plot. 
#' Default to NULL.
#' @param ... further arguments to be passed to the plot() function
#' @return the decoded series of latent states
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' \dontrun{
#' #3 zero-inflated poissons
#' prior_init <- c(0.3,0.3,0.4)
#' dt_init <- c(10,8,6)
#' emit_init <- c(10,50,100)
#' zeroprop <- c(0.5,0.3,0.2)
#' trunc <- c(10,10,10)
#' omega <- matrix(c(0,0.3,0.7,0.4,0,0.6,0.5,0.5,0),3,3,byrow=TRUE)
#' result <- hsmmsim(n=1000,M=3,prior=prior_init,dt_dist="shiftpoisson",
#'          dt_parm=dt_init, tpm_parm=omega,emit_parm=emit_init,zeroprop=zeroprop)
#' y <- result$series
#' state <- result$state
#' fit <- hsmmfit(y=y,ntimes=NULL,M=3,trunc=trunc,prior_init=prior_init,dt_dist="shiftpoisson",
#'      dt_init=dt_init,tpm_init=omega,emit_init=emit_init,zero_init=zeroprop,
#'      method="Nelder-Mead",hessian=FALSE,control=list(maxit=500,trace=1))
#' decode <- hsmmviterbi(y=y,ntimes=NULL,M=3,trunc=trunc,prior=fit$prior,dt_dist="shiftpoisson",
#'      dt_parm=fit$dt_parm,tpm_parm=fit$tpm,emit_parm=fit$emit_parm,
#'      zero_init=fit$zeroprop,plot=TRUE,xlim=c(0,1000),ylim=c(0,200))
#'#check the missclassification rate
#'sum(decode!=state)/length(state)
#'}   
#'      
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hsmmviterbi <- function(y, ntimes=NULL, M, trunc, prior, dt_dist="nonparametric", 
                        dt_parm, tpm_parm, emit_parm, zero_init,
                        plot=TRUE, xlim=NULL, ylim=NULL,...){
  
  if(!dt_dist%in%c("log","shiftpoisson","nonparametric")) 
    stop("dt_dist can only be 'log' or 'shiftpoisson'!")
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(prior)!=M | length(emit_parm)!=M | length(zero_init)!=M |
     nrow(tpm_parm)!= M | ncol(tpm_parm)!=M | length(trunc)!=M) 
    stop("The dimension of the initial value does not equal M!") 
  
    if(is.null(ntimes)) ntimes <- length(y)
    state <- rep(NA, sum(ntimes))
    lastindex <- 0
    for(i in 1:length(ntimes)){
      if(dt_dist=="nonparametric"){
        state[(lastindex+1):(lastindex+ntimes[i])] <- hsmm_viterbi_np(y[(lastindex+1):(lastindex+ntimes[i])],
                                                                      M,trunc,prior,emit_parm,tpm_parm,dt_parm,
                                                                      zero_init)
      }else{
      state[(lastindex+1):(lastindex+ntimes[i])] <- hsmm_viterbi(y[(lastindex+1):(lastindex+ntimes[i])],
                                                                 M, prior, emit_parm, tpm_parm, dt_parm, 
                                                                 trunc, dt_dist, zero_init)
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

