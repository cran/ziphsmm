
########################################################
#' Viterbi algorithm to decode the latent states in hidden semi-Markov models
#' with covariates where the latent state durations have accelerated failure
#' time structure
#' @param y the observed series to be decoded
#' @param M number of latent states
#' @param trunc a vector specifying the truncation at the maximum number of dwelling time 
#' in each state.
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
#' @param plot whether a plot should be returned
#' @param xlim vector specifying the minimum and maximum on the x-axis in the plot. 
#' Default to NULL.
#' @param ylim vector specifying the minimum and maximum on the y-axis in the plot. 
#' Default to NULL.
#' @param ... further arguments to be passed to the plot() function 
#' @return decoded series of latent states
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
#' Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
#' @examples
#' \dontrun{
#' M <- 3
#' prior <- c(0.5,0.3,0.2)
#' dtrate <- c(6,5,4)
#' dtparm <- matrix(c(0.2,0.1,0.2),nrow=3)
#' zeroparm <- c(0,-0.2)
#' emitparm <- matrix(c(4,0.3,5,0.2,6,-0.1),3,2,byrow=TRUE)
#' tpmparm <- c(1,0.2,0.5,-0.2,0,0.2)
#'
#' emit_x <- matrix(c(rep(1,1000),rep(0,1000)),nrow=2000,ncol=1)
#' dt_x <- emit_x
#' tpm_x <- emit_x
#' zeroinfl_x <- emit_x
#' trunc <- c(18,15,10)
#'
#' re <- hsmmsim2_exp(prior,dtrate,dtparm,zeroparm,emitparm,tpmparm,
#'                   trunc, M, n, dt_x,tpm_x, emit_x, zeroinfl_x)
#' y <- re$series
#'
#' rrr <- hsmmfit_exp(y,M,trunc,dtrate,dtparm,prior,zeroparm,emitparm,tpmparm,
#'                   dt_x,zeroinfl_x,emit_x,tpm_x,method="BFGS",control=list(trace=1))
#'
#' decode <- hsmmviterbi_exp(y,M, trunc,dtrate,dtparm,
#'                           prior,zeroparm,emitparm,tpmparm,
#'                           dt_x, zeroinfl_x, emit_x, tpm_x)
#' sum(decode!=re$state)
#' }
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hsmmviterbi_exp <- function(y, M, trunc,dtrate,dtparm,
                            prior,zeroparm,emitparm,tpmparm,
                            dt_x, zeroinfl_x, emit_x, tpm_x,
                            plot=FALSE, xlim=NULL, ylim=NULL, ...){
   
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  n <- length(y)
  ntimes <- n
  state <- rep(NA, n)
  
  tempmat <- matrix(0,nrow=length(y),ncol=1) #for NULL covariates
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
  
  state <- hsmm_cov_viterbi_exp(y,trunc,M,
                                prior,dtrate,dtparm,tpmparm,zeroparm,
                                emitparm,ncovtpm,
                                covdt,covtpm,covzeroinfl,covemit)
  
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
