% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hmmsim2.cont.R
\name{hmmsim2.cont}
\alias{hmmsim2.cont}
\title{Simulate a continuous-time hidden Markov series and its underlying 
states with covariates}
\usage{
hmmsim2.cont(workparm, M, n, zeroindex, emit_x = NULL, zeroinfl_x = NULL,
  timeindex)
}
\arguments{
\item{workparm}{working parameters}

\item{M}{number of latent states}

\item{n}{length of the simulated series}

\item{zeroindex}{a vector specifying whether a certain state is zero-inflated}

\item{emit_x}{matrix of covariates for the log poisson means. Default to NULL.}

\item{zeroinfl_x}{matrix of covariates for the nonzero structural zero proportions. Default to NULL.}

\item{timeindex}{a vector containing the time points}
}
\value{
simulated series and corresponding states
}
\description{
Simulate a continuous-time hidden Markov series and its underlying 
states with covariates
}
\examples{
priorparm <- 0
tpmparm <- c(-1,-2)
zeroindex <- c(1,0)
zeroparm <- c(0,-1,1)
emitparm <- c(2,0.5,-0.5,3,0.3,-0.2)
workparm <- c(priorparm,tpmparm,zeroparm,emitparm)
timeindex <- rep(1,1000)
for(i in 2:1000) timeindex[i] <- timeindex[i-1] + sample(1:4,1)

designx <- matrix(rnorm(2000),nrow=1000,ncol=2)
result <- hmmsim2.cont(workparm,2,1000,zeroindex,emit_x=designx,
                      zeroinfl_x=designx,timeindex=timeindex)
y <- result$series
state <- result$state

}
\references{
Walter Zucchini, Iain L. MacDonald, Roland Langrock. Hidden Markov Models for 
Time Series: An Introduction Using R, Second Edition. Chapman & Hall/CRC
}
