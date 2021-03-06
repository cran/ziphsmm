% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsmmfit.R
\name{hsmmfit}
\alias{hsmmfit}
\title{Estimate the parameters of a general zero-inflated Poisson hidden semi-Markov model
by directly minimizing of the negative log-likelihood function using the gradient
descent algorithm.}
\usage{
hsmmfit(y, ntimes = NULL, M, trunc, prior_init, dt_dist, dt_init, tpm_init,
  emit_init, zero_init, prior_x = NULL, dt_x = NULL, tpm_x = NULL,
  emit_x = NULL, zeroinfl_x = NULL, method = "Nelder-Mead",
  hessian = FALSE, ...)
}
\arguments{
\item{y}{observed time series values}

\item{ntimes}{A vector specifying the lengths of individual, 
i.e. independent, time series. If not specified, the responses are assumed to 
form a single time series, i.e. ntimes=length(data)}

\item{M}{number of hidden states}

\item{trunc}{a vector specifying truncation at the maximum number of dwelling time in each state. 
The higher the truncation, the more accurate the approximation but also the more 
computationally expensive.}

\item{prior_init}{a vector of initial value for prior probability for each state}

\item{dt_dist}{dwell time distribution, can only be "log", "geometric",
or "shiftedpoisson"}

\item{dt_init}{a vector of initial value for the parameter in each dwell time distribution, which 
should be a vector of p's for dt_dist == "log" and a vector of theta's for 
dt_dist=="shiftpoisson"}

\item{tpm_init}{a matrix of initial values for the transition probability matrix, whose diagonal
elements should be zero's}

\item{emit_init}{a vector initial value for the vector containing means for each poisson distribution}

\item{zero_init}{a vector initial value for the vector containing structural zero proportions in each state}

\item{prior_x}{matrix of covariates for generalized logit of prior probabilites (excluding the 
1st probability). Default to NULL.}

\item{dt_x}{matrix of covariates for the dwell time distribution parameters}

\item{tpm_x}{matrix of covariates for transition probability matrix (excluding the 1st column).
Default to NULL.}

\item{emit_x}{matrix of covariates for the log poisson means. Default to NULL.}

\item{zeroinfl_x}{matrix of covariates for the nonzero structural zero proportions. Default to NULL.}

\item{method}{method to be used for direct numeric optimization. See details in
the help page for optim() function. Default to Nelder-Mead.}

\item{hessian}{Logical. Should a numerically differentiated Hessian matrix be returned?
Note that the hessian is for the working parameters, which are the logit of parameter p for 
each log-series dwell time distribution or the log of parameter theta for each 
shifted-poisson dwell time distribution, the generalized logit of prior probabilities (except for the
1st state),the logit of each nonzero structural zero proportions, the log of each 
state-dependent poisson means, and the generalized logit of the transition probability 
matrix(except 1st column and the diagonal elements)}

\item{...}{Further arguments passed on to the optimization methods}
}
\value{
simulated series and corresponding states
}
\description{
Estimate the parameters of a general zero-inflated Poisson hidden semi-Markov model
by directly minimizing of the negative log-likelihood function using the gradient
descent algorithm.
}
\examples{

#2 zero-inflated poissons
prior_init <- c(0.5,0.5)
emit_init <- c(10,30)
dt_init <- c(10,6)
trunc <- c(20,10)
zeroprop <- c(0.5,0.3)
omega <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
sim2 <- hsmmsim(n=1000,M=2,prior=prior_init,dt_dist="shiftpoisson",
         dt_parm=dt_init, tpm_parm=omega,
         emit_parm=emit_init,zeroprop=zeroprop)
str(sim2)
y <- sim2$series
fit2 <- hsmmfit(y=y,M=2,trunc=trunc,prior_init=prior_init,dt_dist="shiftpoisson",
     dt_init=dt_init,
     tpm_init=omega,emit_init=emit_init,zero_init=zeroprop,
     method="Nelder-Mead",hessian=FALSE,control=list(maxit=500,trace=1))
str(fit2)


\dontrun{
#1 zero-inflated poisson and 3 regular poissons
prior_init <- c(0.5,0.2,0.2,0.1)
dt_init <- c(0.8,0.7,0.6,0.5)
emit_init <- c(10,30,70,130)
trunc <- c(10,10,10,10)
zeroprop <- c(0.6,0,0,0)  #only the 1st-state is zero-inflated
omega <- matrix(c(0,0.5,0.3,0.2,0.4,0,0.4,0.2,
                  0.2,0.6,0,0.2,0.1,0.1,0.8,0),4,4,byrow=TRUE)
sim1 <- hsmmsim(n=2000,M=4,prior=prior_init,dt_dist="log",
         dt_parm=dt_init, tpm_parm=omega,
         emit_parm=emit_init,zeroprop=zeroprop)
str(sim1)
y <- sim1$series
fit <- hsmmfit(y=y,M=4,trunc=trunc,prior_init=prior_init,dt_dist="log",dt_init=dt_init,
     tpm_init=omega,emit_init=emit_init,zero_init=zeroprop,
     method="Nelder-Mead",hessian=TRUE,control=list(maxit=500,trace=1))
str(fit)

#variances for the 20 working parameters, which are the logit of parameter p for 
#the 4 log-series dwell time distributions, the generalized logit of prior probabilities 
#for state 2,3,4, the logit of each nonzero structural zero proportions in state 1, 
#the log of 4 state-dependent poisson means, and the generalized logit of the 
#transition probability matrix(which are tpm[1,3],tpm[1,4], tpm[2,3],tpm[2,4],
#tpm[3,2],tpm[3,4],tpm[4,2],tpm[4,3])
variance <- diag(solve(fit$obsinfo)) 


#1 zero-inflated poisson and 2 poissons with covariates
data(CAT)
y <- CAT$activity
x <- data.matrix(CAT$night)
prior_init <- c(0.5,0.3,0.2)
dt_init <- c(0.9,0.6,0.3)
emit_init <- c(10,20,30)
zero_init <- c(0.5,0,0) #assuming only the 1st state has structural zero's
tpm_init <- matrix(c(0,0.3,0.7,0.4,0,0.6,0.5,0.5,0),3,3,byrow=TRUE)
trunc <- c(10,7,4)
fit2 <-  hsmmfit(y,rep(1440,3),3,trunc,prior_init,"log",dt_init,tpm_init,
     emit_init,zero_init,emit_x=x,zeroinfl_x=x,hessian=FALSE,
     method="Nelder-Mead", control=list(maxit=500,trace=1))
fit2

#another example with covariates for 2 zero-inflated poissons
data(CAT)
y <- CAT$activity
x <- data.matrix(CAT$night)
prior_init <- c(0.5,0.5)
dt_init <- c(10,5)
emit_init <- c(10, 30)
zero_init <- c(0.5,0.2)
tpm_init <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
trunc <- c(10,5)
fit <-  hsmmfit(y,NULL,2,trunc,prior_init,"shiftpoisson",dt_init,tpm_init,
     emit_init,zero_init,dt_x=x,emit_x=x,zeroinfl_x=x,tpm_x=x,hessian=FALSE,
     method="Nelder-Mead", control=list(maxit=500,trace=1))
fit
}

}
\references{
Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
}
