
########################################################
#' Viterbi algorithm to decode the latent states in continuous-time
#' hidden Markov models with covariates
#' @param y the observed series to be decoded
#' @param M number of latent states
#' @param workparm a vector of values for working parameters, which is the last element returned 
#' from hmmfit() function. This consists the generalized logit of prior probabilities 
#' (except for the 1st state), generalized logit of transition probability matrix
#' (except for the 1st column), the  logit of nonzero structural zero proportions, and 
#' the log poisson means
#' @param zero_init a vector containing structural zero proportions in each state, e.g. set
#' zero_init[i] to be 0 if the i-th state is a regular poisson, and otherwise 1.
#' @param emit_x matrix of covariates for the log poisson means. Default to NULL.
#' @param zeroinfl_x matrix of covariates for the nonzero structural zero proportions. Default to NULL.
#' @param timeindex a vector containing the time points
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
#' priorparm <- 0
#' tpmparm <- c(-1,-2)
#' zeroindex <- c(1,0)
#' zeroparm <- c(0,-1,1)
#' emitparm <- c(2,0.5,-0.5,3,0.3,-0.2)
#' workparm <- c(priorparm,tpmparm,zeroparm,emitparm)
#' timeindex <- rep(1,1000)
#' for(i in 2:1000) timeindex[i] <- timeindex[i-1] + sample(1:4,1)
#'
#' designx <- matrix(rnorm(2000),nrow=1000,ncol=2)
#' x <- cbind(1,designx) #has to make the additional 1st column of 1 for intercept
#' result <- hmmsim2.cont(workparm,2,1000,zeroindex,emit_x=designx,
#'                       zeroinfl_x=designx,timeindex=timeindex)
#' y <- result$series
#' state <- result$state
#'
#' fit2 <-  fasthmmfit.cont(y=y,x=designx,M=2,prior_init=c(0.5,0.5),
#'   tpm_init=matrix(c(-0.2,0.2,0.1,-0.1),2,2,byrow=TRUE),
#'   zero_init=0.4,emit_init=c(7,21), timeindex=timeindex,
#'   hessian=FALSE, method="BFGS", control=list(trace=1))

#' decode2 <- hmmviterbi2.cont(y,2,fit2$working_parm,c(1,0),
#'   emit_x=designx, zeroinfl_x=designx,
#'   timeindex=timeindex,plot=FALSE)
#' @useDynLib ziphsmm
#' @importFrom Rcpp evalCpp
#' @export
hmmviterbi2.cont <- function(y, M, workparm, zero_init,
                        emit_x=NULL,zeroinfl_x=NULL,timeindex,
                        plot=FALSE, xlim=NULL, ylim=NULL, ...){
  
  if(floor(M)!=M | M<2) stop("The number of latent states must be an integer greater than or equal to 2!")
  if(length(zero_init)!=M) stop("The dimension of zero_init is not M!")
  ntimes <- length(y)
  
  state <- rep(NA, sum(ntimes))
  
  tempmat <- matrix(0,nrow=length(y),ncol=1) #for NULL covariates
  
    ncovprior <- 0
    covprior <- tempmat
 
    ncovtpm <- 0
    covtpm <- tempmat

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

    vdiff <- diff(timeindex)
    udiff <- sort(unique(vdiff))
    tpm_init <- retrieve_cov_cont(workparm,M,ncovemit+1)$gamma
    expms <- getallexpm(tpm_init, udiff)
    
  state <- hmm_cov_viterbi_cont(workparm, M,y,ncovprior,covprior,
                          ncovtpm,covtpm,ncovzeroinfl,covzeroinfl,
                          ncovemit,covemit, zero_init, timeindex, udiff, expms)

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
