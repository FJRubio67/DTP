% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DTP.R
\name{ddtp}
\alias{ddtp}
\title{Probability Density Function for DTP distributions}
\usage{
ddtp(x, mu, par1, par2, delta1, delta2, f, param = "tp", log = FALSE)
}
\arguments{
\item{log, }{log.p: logical; if TRUE, probabilities p are given as log(p).}

\item{x:}{vector of quantiles.}

\item{p:}{vector of probabilities.}

\item{n:}{number of observations. If length(n) > 1, the length is taken to be the number required.}

\item{mu:}{location parameter.}

\item{par1:}{scale parameter 1.}

\item{par2:}{scale parameter 2.}

\item{delta1:}{shape parameter 1.}

\item{delta2:}{shape parameter 2.}

\item{F, }{qF, rF, f: distribution function, quantile function, random function and density function associated of a symmetric random variable.}

\item{param:}{parameterisations used.}
}
\value{

}
\description{
Probability Density Function for DTP distributions
}
