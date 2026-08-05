rm(list=ls())


#********************************************************************************
# Baseline functions
#********************************************************************************


#######################################################################
# Symmetric sinh-arcsinh (SAS) distribution
#######################################################################

# baseline symmetric SAS pdf
dsas0 <- function(x,delta,log=FALSE){
  logPDF <- dnorm(sinh(delta*asinh(x)),log=T) + log(delta) + log(cosh(delta*asinh(x))) -0.5*log(1+x^2)
  ifelse( is.numeric(logPDF),ifelse( log, return(logPDF), return(exp(logPDF)) ), logPDF )
}

# baseline symmetric SAS cdf
psas0 <- function(x,delta,log.p=FALSE){
  logCDF <- pnorm(sinh(delta*asinh(x)),log.p=T)
  ifelse( is.numeric(logCDF),ifelse( log.p, return(logCDF), return(exp(logCDF)) ), logCDF )
}


#####################################################
# Other functions
#####################################################

# logit function
logit <- Vectorize(function(p){
  val <- log(p) -  log(1-p)
  return(as.vector(val))
})

# expit functon
expit <- Vectorize(function(x){
  a <- exp(x)
  return(as.vector( a/(a + 1) ))          
})





# Required packages
library(mnormt)
library(mvtnorm)
library(ghyp)
library(copula)
library(spBayes)
library(fMultivar)
# library(devtools)
# github_install("FJRubio67/twopiece")
library(twopiece)
# github_install("FJRubio67/DTP")
library(DTP)

# routines
source("/Users/FJRubio/Dropbox/Ia_Ro_Ru/code/routines.R")


# The SMI stocks data:  "SMI"  vs    "Swiss.Re"
data(smi.stocks)
Y = as.matrix(smi.stocks[ , c(1, 6)])
colnames(smi.stocks[ , c(1, 6)])

plot(Y)

##############################################################################################
# bivariate Gaussian-copula with DTP SAS marginals
##############################################################################################

# log likelihood function (reparameterised)
loglikDTP = function(par){
  # Reparameterisation
  mu = par[1:2]; sigma = exp(par[3:4]); 
  epsilon = 2*expit(par[5:6])-1; delta1 = exp(par[7:8]); delta2 <- exp(par[9:10])
  rho = 2*expit(par[11])-1; 
  # Ingredients of the log likelihood
  norm.cop = normalCopula(rho, dim = 2)
  probs =cbind(pdtp(Y[,1],mu[1],sigma[1],epsilon[1],delta1[1],delta1[2],param="eps",F=psas0, f=dsas0),
               pdtp(Y[,2],mu[2],sigma[2],epsilon[2],delta2[1],delta2[2],param="eps",F=psas0, f=dsas0 ))
  val.cop = sum(dCopula(probs, copula = norm.cop, log = TRUE))
  val.marg1 = sum( ddtp(Y[,1],mu[1],sigma[1],epsilon[1],delta1[1],delta1[2],param="eps",f=dsas0,log=T) )
  val.marg2 = sum( ddtp(Y[,2],mu[2],sigma[2],epsilon[2],delta2[1],delta2[2],param="eps",f=dsas0,log=T) )
  
  # output
  out <- -val.cop - val.marg1 - val.marg2
  
  return(out)
}


start2 = c(0,0,log(0.01),log(0.01),0,0,-0.5,-0.5,-0.5,-0.5,0) 

# Optimisation step
OPT.DTP = optim(start2,loglikDTP,control=list(maxit=10000))
OPT.DTP

MLE.DTP <- c(OPT.DTP$par[1:2], exp(OPT.DTP$par[3:4]), 2*expit(OPT.DTP$par[5:6])-1, exp(OPT.DTP$par[7:10]),  2*expit(OPT.DTP$par[11])-1  )
MLE.DTP
AIC.DTP <- 2*OPT.DTP$value + 2*length(MLE.DTP)
BIC.DTP <- 2*OPT.DTP$value + log(length(nrow(Y)))*length(MLE.DTP)

##############################################################################################
# bivariate Gaussian-copula with TPSC SAS marginals
##############################################################################################

# log likelihood function (reparameterised)
loglikTPSC = function(par){
  # Reparameterisation
  mu = par[1:2]; sigma = exp(par[3:4]); 
  epsilon = 2*expit(par[5:6])-1; delta = exp(par[7:8]);
  rho = 2*expit(par[9])-1; 
# Ingredients of the log likelihood
    norm.cop = normalCopula(rho, dim = 2)
    probs =cbind(ptp4(Y[,1],mu[1],sigma[1],epsilon[1],delta[1],param="eps",FUN=psas0),
                 ptp4(Y[,2],mu[2],sigma[2],epsilon[2],delta[2],param="eps",FUN=psas0))
    val.cop = sum(dCopula(probs, copula = norm.cop, log = TRUE))
    val.marg1 = sum( dtp4(Y[,1],mu[1],sigma[1],epsilon[1],delta[1],param="eps",FUN=dsas0,log=T) )
    val.marg2 = sum( dtp4(Y[,2],mu[2],sigma[2],epsilon[2],delta[2],param="eps",FUN=dsas0,log=T) )

    # output
    out <- -val.cop - val.marg1 - val.marg2
    
        return(out)
}


start = c(0,0,log(0.01),log(0.01),0,0,-0.5,-0.5,0) 

# Optimisation step
OPT.TPSC = optim(start,loglikTPSC,control=list(maxit=10000))
OPT.TPSC

MLE.TPSC <- c(OPT.TPSC$par[1:2], exp(OPT.TPSC$par[3:4]), 2*expit(OPT.TPSC$par[5:6])-1, exp(OPT.TPSC$par[7:8]),  2*expit(OPT.TPSC$par[9])-1  )
MLE.TPSC
AIC.TPSC <- 2*OPT.TPSC$value + 2*length(MLE.TPSC)
BIC.TPSC <- 2*OPT.TPSC$value + log(length(nrow(Y)))*length(MLE.TPSC)
