#############################################################################################################3
# mu   : location parameter 
# sigma1   : scale parameter 1
# sigma1   : scale parameter 1
# delta1      : shape parameter 1
# delta2    : shape parameter 2 
#############################################################################################################3
rm(list=ls())
library(DTP)
#********************************************************************************
# Baseline functions
#********************************************************************************
# baseline SAS pdf
dsas0 <- function(x,delta,log=FALSE){
  logPDF <- dnorm(sinh(delta*asinh(x)),log=T) + log(delta) + log(cosh(delta*asinh(x))) -0.5*log(1+x^2)
  ifelse( is.numeric(logPDF),ifelse( log, return(logPDF), return(exp(logPDF)) ), logPDF )
}

# baseline SAS cdf
psas0 <- function(x,delta,log.p=FALSE){
  logCDF <- pnorm(sinh(delta*asinh(x)),log.p=T)
  ifelse( is.numeric(logCDF),ifelse( log.p, return(logCDF), return(exp(logCDF)) ), logCDF )
}

# baseline quantile function
qsas0 <- function(p,delta){
  Q <- sinh((asinh(qnorm(p)))/delta)
  return(Q)
}

# baseline RNG
rsas0 <- function(n,delta){
  sample <- sinh((asinh(rnorm(n)))/delta)
  return(sample)
}

#********************************************************************************
# DTP functions
#********************************************************************************
# Probability Density Function
ddtpsas <- function(x, mu, sigma1, sigma2, delta1, delta2, log = FALSE){
  out <- ddtp(x, mu, sigma1, sigma2, delta1, delta2, dsas0, param = "tp", log = log)
  return(out)
}

# Cumulative Distribution Function
pdtpsas <- function(x, mu, sigma1, sigma2, delta1, delta2, log.p = FALSE){
  out <- pdtp(x, mu, sigma1, sigma2, delta1, delta2, psas0, dsas0, param = "tp", log = FALSE)
  return(out)
}

# Quantile Function
qdtpsas <- function(x, mu, sigma1, sigma2, delta1, delta2, log.p = FALSE){
  out <- qdtp(x, mu, sigma1, sigma2, delta1, delta2, rsas0, dsas0, param = "tp", log = FALSE)
  return(out)
}

# RNG Function
rdtpsas <- function(x, mu, sigma1, sigma2, delta1, delta2, log.p = FALSE){
  out <- rdtp(x, mu, sigma1, sigma2, delta1, delta2, rsas0, dsas0, param = "tp")
  return(out)
}




#############################################################################################################3
# Illustrations
#############################################################################################################3

# Simulated data
set.seed(123)
data <- rdtpsas(1000,0,1,2,0.75,1.5)
# True pdf
pdf.true <- Vectorize(function(x) ddtpsas(x,0,1,2,0.75,1.5))
# True cdf
cdf.true <- Vectorize(function(x) pdtpsas(x,0,1,2,0.75,1.5))


# Histogram vs pdf
hist(data, breaks = 30, probability = T, cex.axis = 1.5, cex.lab = 1.5, 
     xlab = "x", ylab = "density", main = "Histogram vs pdf")
curve(pdf.true, -5,15, add= T, col="red", lwd = 2, n = 1000)
box()


# ECDF vs cdf
plot(ecdf(data), cex.axis = 1.5, cex.lab = 1.5, 
     xlab = "x", ylab = "density", main = "ECDF vs cdf")
curve(cdf.true, -5,15, add= T, col="gray", lwd = 2, n = 1000)
box()