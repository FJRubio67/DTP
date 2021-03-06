% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DTP.R
\name{qdtp}
\alias{qdtp}
\title{Quantile Function for DTP distributions}
\usage{
qdtp(p, mu, par1, par2, delta1, delta2, qF, f, param = "tp")
}
\arguments{
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

\item{log, }{log.p: logical; if TRUE, probabilities p are given as log(p).}
}
\value{

}
\description{
Quantile Function for DTP distributions
}
