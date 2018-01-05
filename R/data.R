#' Pseudo activity counts (per minute) data for cats
#'
#' @format A data frame with 4320 rows and 5 variables:
#' \describe{
#'   \item{id}{cat ID: 1,2,3}
#'   \item{hour}{hour of the day: 1,2,...,24}
#'   \item{minute}{minute of the hour: 1,2,...,60}
#'   \item{night}{night time indicator}
#'   \item{activity}{activity count data}
#' }
"CAT"


#prior_init <- c(0.33,0.33,0.33)
#emit_init <- c(20,80,150)
#zeroprop <- c(0.3,0,0)
#omega <- matrix(c(0.5,0.3,0.2,0.4,0.3,0.3,0.2,0.4,0.4),3,3,byrow=TRUE)
#result1 <- hmmsim(n=360,M=3,prior=prior_init, tpm_parm=omega,
#           emit_parm=emit_init,zeroprop=zeroprop)
#emit_init <- c(10,50,100)
#zeroprop <- c(0.6,0,0)
#result2 <- hmmsim(n=720,M=3,prior=prior_init, tpm_parm=omega,
#                   emit_parm=emit_init,zeroprop=zeroprop)
#emit_init <- c(20,80,150)
#zeroprop <- c(0.3,0,0)
#result3 <- hmmsim(n=720,M=3,prior=prior_init, tpm_parm=omega,
#                   emit_parm=emit_init,zeroprop=zeroprop)
#emit_init <- c(10,50,100)
#zeroprop <- c(0.6,0,0)
#result4 <- hmmsim(n=720,M=3,prior=prior_init, tpm_parm=omega,
#                    emit_parm=emit_init,zeroprop=zeroprop)
#emit_init <- c(20,80,150)
#zeroprop <- c(0.3,0,0)
#result5 <- hmmsim(n=720,M=3,prior=prior_init, tpm_parm=omega,
#                   emit_parm=emit_init,zeroprop=zeroprop)
#emit_init <- c(10,50,100)
#zeroprop <- c(0.6,0,0)
#result6 <- hmmsim(n=720,M=3,prior=prior_init, tpm_parm=omega,
#                   emit_parm=emit_init,zeroprop=zeroprop)
#emit_init <- c(20,80,150)
#zeroprop <- c(0.3,0,0)
#omega <- matrix(c(0.5,0.3,0.2,0.4,0.3,0.3,0.2,0.4,0.4),3,3,byrow=TRUE)
#result7 <- hmmsim(n=360,M=3,prior=prior_init, tpm_parm=omega,
#                   emit_parm=emit_init,zeroprop=zeroprop)


#id <- rep(c(1,2,3),c(1440,1440,1440))
#hour <- merge$hour[1:4320]
#minute <- merge$minute[1:4320]
#night <- merge$night[1:4320]
#activity <- c(result1$series,result2$series,result3$series,result4$series,
#              result5$series,result6$series,result7$series)
#CAT <- data.frame(cbind(id,hour,minute,night,activity))
#save(CAT,file="~/Desktop/ziphsmm/data/CAT.RData")
