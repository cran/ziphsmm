
########################################################
#' Viterbi algorithm to decode the latent states for hidden Markov models
#' @param y the observed series to be decoded
#' @param ntimes vector specifying the lengths of individual, 
#' i.e. independent, time series. If not specified, the responses are assumed to 
#' form a single time series, i.e. ntimes=length(y)
#' @param M number of latent states
#' @param prior_init a vector of prior probability values
#' @param tpm_init transition probability matrix
#' @param emit_init a vector containing means for each poisson distribution
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
#' prior_init <- c(0.5,0.2,0.3)
#' emit_init <- c(10,40,70)
#' zero_init <- c(0.5,0,0)
#' omega <- matrix(c(0.5,0.3,0.2,0.4,0.3,0.3,0.2,0.4,0.4),3,3,byrow=TRUE)
#' result <- hmmsim(n=1000,M=3,prior=prior_init, tpm_parm=omega,
#'          emit_parm=emit_init,zeroprop=zero_init)
#' y <- result$series
#' state <- result$state
#' fit <- hmmfit(y=y,M=3,prior_init=prior_init,tpm_init=omega,
#'      emit_init=emit_init,zero_init=zero_init,
#'      method="Nelder-Mead",hessian=FALSE,control=list(maxit=500,trace=1))
#' 
#' decode <- hmmviterbi(y,NULL,3,fit$prior,fit$tpm,fit$emit_parm,fit$zeroprop,
#'                       plot=TRUE,xlab="Time",ylab="Count")      
#' #check the missclassification rate
#' sum(decode!=state)/length(state)      
#' 
#' \dontrun{
#' decode <- hmmviterbi(y,NULL,3,fit$prior,fit$tpm,fit$emit_parm,fit$zeroprop,
#'                       plot=TRUE,xlim=c(0,100),ylim=c(0,100),
#'                       xlab="Time",ylab="Count")  
#' }
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hmmviterbi <- function(y, ntimes=NULL, M, prior_init, tpm_init, 
                          emit_init, zero_init, plot=FALSE,
                       xlim=NULL, ylim=NULL, ...){
 
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(prior_init)!=M | length(emit_init)!=M | length(zero_init)!=M |
     nrow(tpm_init)!= M | ncol(tpm_init)!=M) stop("The dimension of the initial value does not equal M!") 
  
    if(is.null(ntimes)) ntimes <- length(y)
    state <- rep(NA, sum(ntimes))
    lastindex <- 0
    for(i in 1:length(ntimes)){
      state[(lastindex+1):(lastindex+ntimes[i])] <- hmm_viterbi(prior_init,tpm_init,
                                                        emit_init,M,
                                                        y[(lastindex+1):(lastindex+ntimes[i])],
                                                        zero_init)
      
      lastindex <- lastindex+ntimes[i]
    }
    
    print("Be careful of label switching when interpreting states!")
    #possible label switching
    
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