
########################################################
#' Viterbi algorithm to decode the latent states for continuous-time
#' hidden Markov models without covariates
#' @param y the observed series to be decoded
#' @param M number of latent states
#' @param prior_init a vector of prior probability values
#' @param tpm_init transition rate matrix
#' @param emit_init a vector containing means for each poisson distribution
#' @param zero_init a vector containing structural zero proportions in each state
#' @param timeindex a vector containing the time points
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
#' omega <- matrix(c(-0.3,0.2,0.1,0.1,-0.2,0.1,0.2,0.2,-0.4),3,3,byrow=TRUE)
#' timeindex <- rep(1,1000)
#' for(i in 2:1000) timeindex[i] <- timeindex[i-1] + sample(1:3,1)
#' result <- hmmsim.cont(n=1000,M=3,prior=prior_init, tpm_parm=omega,
#'           emit_parm=emit_init,zeroprop=zero_init,timeindex=timeindex)
#' y <- result$series
#' fit2 <-  fasthmmfit.cont(y,x=NULL,M=3,prior_init,omega,
#'                         emit_init,0.5,timeindex=timeindex,hessian=FALSE,
#'                         method="BFGS", control=list(maxit=500,trace=1))
#' decode2 <- hmmviterbi.cont(y,3,fit2$prior,fit2$tpm,fit2$emit,
#'            c(fit2$zeroprop,0,0),timeindex=timeindex)    
#'
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hmmviterbi.cont <- function(y, M, prior_init, tpm_init, 
                       emit_init, zero_init, timeindex, plot=FALSE,
                       xlim=NULL, ylim=NULL, ...){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(prior_init)!=M | length(emit_init)!=M | length(zero_init)!=M |
     nrow(tpm_init)!= M | ncol(tpm_init)!=M) stop("The dimension of the initial value does not equal M!") 
  
  ntimes <- length(y)
  vdiff <- diff(timeindex)
  udiff <- sort(unique(vdiff))
  expms <- getallexpm(tpm_init, udiff)
  
  state <- rep(NA, ntimes)
  lastindex <- 0 
  state <- hmm_viterbi_cont(prior_init,tpm_init,emit_init,M,
                       y,zero_init,timeindex,udiff,expms)
  
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