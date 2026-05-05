#' Probability Density Function for DTP distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param delta1: shape parameter 1.
#' @param delta2: shape parameter 2.
#' @param F, qF, rF, f: distribution function, quantile function, random function and density function associated of a symmetric random variable.
#' @param param: parameterisations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
ddtp <- function (x, mu, par1, par2, delta1, delta2, f, param = "tp",
                  log = FALSE)
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  # FIX (item 2): replaced ifelse()-based validation with proper if/stop guards.
  # FIX (item 2): removed trailing ifelse(is.numeric(logPDF), ...) workaround.
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 or/and par2 or/and delta1 or/and delta2 out of range in the parametrization tp")
    }
    sigma1 = par1
    sigma2 = par2
    logeps   <- log(sigma1) + f(0, delta2, log = TRUE) - log(sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    log1eps  <- log(sigma2) + f(0, delta1, log = TRUE) - log(sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    logPDF <- log(2) + ifelse(x < mu,
                              logeps  + f((x - mu)/sigma1, delta1, log = TRUE) - log(sigma1),
                              log1eps + f((x - mu)/sigma2, delta2, log = TRUE) - log(sigma2))
  }
  if (param == "eps") {
    if (!(par1 > 0 & abs(par2) < 1 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 or/and delta1 or/and delta2 is/are not positive or/and abs(par2) out of range in the parametrization eps")
    }
    sigma1 = par1 * (1 + par2)
    sigma2 = par1 * (1 - par2)
    logeps   <- log(sigma1) + f(0, delta2, log = TRUE) - log(sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    log1eps  <- log(sigma2) + f(0, delta1, log = TRUE) - log(sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    logPDF <- log(2) + ifelse(x < mu,
                              logeps  + f((x - mu)/sigma1, delta1, log = TRUE) - log(sigma1),
                              log1eps + f((x - mu)/sigma2, delta2, log = TRUE) - log(sigma2))
  }
  if (param == "isf") {
    if (!(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 or/and par2 or/and delta1 or/and delta2 out of range in the parametrization isf")
    }
    sigma1 = par1 * par2
    sigma2 = par1/par2
    logeps   <- log(sigma1) + f(0, delta2, log = TRUE) - log(sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    log1eps  <- log(sigma2) + f(0, delta1, log = TRUE) - log(sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    logPDF <- log(2) + ifelse(x < mu,
                              logeps  + f((x - mu)/sigma1, delta1, log = TRUE) - log(sigma1),
                              log1eps + f((x - mu)/sigma2, delta2, log = TRUE) - log(sigma2))
  }
  if (log) return(logPDF) else return(exp(logPDF))
}


#' Cumulative Probability Function for DTP distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param delta1: shape parameter 1.
#' @param delta2: shape parameter 2.
#' @param F, qF, rF, f: distribution function, quantile function, random function and density function associated of a symmetric random variable.
#' @param param: parameterisations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
pdtp <- function (x, mu, par1, par2, delta1, delta2, F, f, param = "tp",
                  log.p = FALSE)
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  # FIX (item 2): replaced ifelse()-based validation with proper if/stop guards.
  # FIX (item 2): removed trailing ifelse(is.numeric(CDF), ...) workaround.
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp")
    }
    sigma1 = par1
    sigma2 = par2
    eps <- sigma1 * f(0, delta2, log = FALSE) / (sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    CDF <- ifelse(x < mu,
                  2 * eps * F((x - mu)/sigma1, delta1, log.p = FALSE),
                  eps + (1 - eps) * (2 * F((x - mu)/sigma2, delta2, log.p = FALSE) - 1))
  }
  if (param == "eps") {
    if (!(par1 > 0 & abs(par2) < 1 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    sigma1 = par1 * (1 + par2)
    sigma2 = par1 * (1 - par2)
    eps <- sigma1 * f(0, delta2, log = FALSE) / (sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    CDF <- ifelse(x < mu,
                  2 * eps * F((x - mu)/sigma1, delta1, log.p = FALSE),
                  eps + (1 - eps) * (2 * F((x - mu)/sigma2, delta2, log.p = FALSE) - 1))
  }
  if (param == "isf") {
    if (!(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    sigma1 = par1 * par2
    sigma2 = par1/par2
    eps <- sigma1 * f(0, delta2, log = FALSE) / (sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    CDF <- ifelse(x < mu,
                  2 * eps * F((x - mu)/sigma1, delta1, log.p = FALSE),
                  eps + (1 - eps) * (2 * F((x - mu)/sigma2, delta2, log.p = FALSE) - 1))
  }
  if (log.p) return(log(CDF)) else return(CDF)
}


#' Random Number Generation for DTP distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param delta1: shape parameter 1.
#' @param delta2: shape parameter 2.
#' @param F, qF, rF, f: distribution function, quantile function, random function and density function associated of a symmetric random variable.
#' @param param: parameterisations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
rdtp <- function (n, mu, par1, par2, delta1, delta2, rF, f, param = "tp")
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  # FIX (item 2): replaced ifelse()-based validation with proper if/stop guards.
  # FIX (item 3): pre-compute a single shared pool of n draws from rF for each
  #               arm so ifelse does not trigger two independent calls to rF(n, ...).
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp")
    }
    sigma1 = par1
    sigma2 = par2
    eps <- sigma1 * f(0, delta2, log = FALSE) / (sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    u  <- runif(n)
    z1 <- abs(rF(n, delta1))
    z2 <- abs(rF(n, delta2))
    draws <- ifelse(u < eps, mu - sigma1 * z1, mu + sigma2 * z2)
  }
  if (param == "eps") {
    if (!(par1 > 0 & abs(par2) < 1 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    sigma1 = par1 * (1 + par2)
    sigma2 = par1 * (1 - par2)
    eps <- sigma1 * f(0, delta2, log = FALSE) / (sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    u  <- runif(n)
    z1 <- abs(rF(n, delta1))
    z2 <- abs(rF(n, delta2))
    draws <- ifelse(u < eps, mu - sigma1 * z1, mu + sigma2 * z2)
  }
  if (param == "isf") {
    if (!(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    sigma1 = par1 * par2
    sigma2 = par1/par2
    eps <- sigma1 * f(0, delta2, log = FALSE) / (sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    u  <- runif(n)
    z1 <- abs(rF(n, delta1))
    z2 <- abs(rF(n, delta2))
    draws <- ifelse(u < eps, mu - sigma1 * z1, mu + sigma2 * z2)
  }
  return(draws)
}



#' Quantile Function for DTP distributions
#' @param x: vector of quantiles.
#' @param p: vector of probabilities.
#' @param n: number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu: location parameter.
#' @param par1: scale parameter 1.
#' @param par2: scale parameter 2.
#' @param delta1: shape parameter 1.
#' @param delta2: shape parameter 2.
#' @param F, qF, rF, f: distribution function, quantile function, random function and density function associated of a symmetric random variable.
#' @param param: parameterisations used.
#' @param log, log.p: logical; if TRUE, probabilities p are given as log(p).
#' @return
#' @export
qdtp <- function (p, mu, par1, par2, delta1, delta2, qF, f, param = "tp")
{
  param = match.arg(param, choices = c("tp", "eps", "isf"))
  # FIX (item 2): replaced ifelse()-based validation with proper if/stop guards.
  if (param == "tp") {
    if (!(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp")
    }
    sigma1 = par1
    sigma2 = par2
    eps <- sigma1 * f(0, delta2, log = FALSE) / (sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    Q <- ifelse(p < eps,
                mu + sigma1 * qF(p/(2 * eps), delta1),
                mu + sigma2 * qF(0.5 * (1 + (p - eps)/(1 - eps)), delta2))
  }
  if (param == "eps") {
    if (!(par1 > 0 & abs(par2) < 1 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
    }
    sigma1 = par1 * (1 + par2)
    sigma2 = par1 * (1 - par2)
    eps <- sigma1 * f(0, delta2, log = FALSE) / (sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    Q <- ifelse(p < eps,
                mu + sigma1 * qF(p/(2 * eps), delta1),
                mu + sigma2 * qF(0.5 * (1 + (p - eps)/(1 - eps)), delta2))
  }
  if (param == "isf") {
    if (!(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0)) {
      stop("invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
    }
    sigma1 = par1 * par2
    sigma2 = par1/par2
    eps <- sigma1 * f(0, delta2, log = FALSE) / (sigma1 * f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, log = FALSE))
    Q <- ifelse(p < eps,
                mu + sigma1 * qF(p/(2 * eps), delta1),
                mu + sigma2 * qF(0.5 * (1 + (p - eps)/(1 - eps)), delta2))
  }
  return(Q)
}
