# Probability Density Function

ddtp <- function (x, mu, par1, par2, delta1, delta2, f, param = "tp", 
                  log = FALSE) 
{
  if (param == "tp") {
    sigma1 = par1
    sigma2 = par2
    logeps <- log(sigma1) + f(0, delta2, log = T) - log(sigma1 * 
                                                          f(0, delta2, log = F) + sigma2 * f(0, delta1, log = F))
    log1eps <- log(sigma2) + f(0, delta1, log = T) - log(sigma1 * 
                                                           f(0, delta2, log = F) + sigma2 * f(0, delta1, log = F))
    ifelse(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0, 
           logPDF <- log(2) + ifelse(x < mu, logeps + f((x - 
                                                           mu)/sigma1, delta1, log = T) - log(sigma1), 
                                     log1eps + f((x - mu)/sigma2, delta2, log = T) - 
                                       log(sigma2)), logPDF <- "invalid arguments: par1 or/and par2 or/and delta1 or/and delta2 out of range in the parametrization tp")
  }
  if (param == "eps") {
    sigma1 = par1 * (1 + par2)
    sigma2 = par1 * (1 - par2)
    logeps <- log(sigma1) + f(0, delta2, log = T) - log(sigma1 * 
                                                          f(0, delta2, log = F) + sigma2 * f(0, delta1, log = F))
    log1eps <- log(sigma2) + f(0, delta1, log = T) - log(sigma1 * 
                                                           f(0, delta2, log = F) + sigma2 * f(0, delta1, log = F))
    ifelse(par1 > 0 & abs(par2) < 1 & delta1 > 0 & delta2 > 
             0, logPDF <- log(2) + ifelse(x < mu, logeps + f((x - 
                                                                mu)/sigma1, delta1, log = T) - log(sigma1), log1eps + 
                                            f((x - mu)/sigma2, delta2, log = T) - log(sigma2)), 
           logPDF <- "invalid arguments: par1 or/and delta1 or/and delta2 is/are not positive or/and abs(par2) out of range in the parametrization eps")
  }
  if (param == "isf") {
    sigma1 = par1 * par2
    sigma2 = par1/par2
    logeps <- log(sigma1) + f(0, delta2, log = T) - log(sigma1 * 
                                                          f(0, delta2, log = F) + sigma2 * f(0, delta1, log = F))
    log1eps <- log(sigma2) + f(0, delta1, log = T) - log(sigma1 * 
                                                           f(0, delta2, log = F) + sigma2 * f(0, delta1, log = F))
    ifelse(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0, 
           logPDF <- log(2) + ifelse(x < mu, logeps + f((x - 
                                                           mu)/sigma1, delta1, log = T) - log(sigma1), 
                                     log1eps + f((x - mu)/sigma2, delta2, log = T) - 
                                       log(sigma2)), logPDF <- "invalid arguments: par1 or/and par2 or/and delta1 or/and delta2 out of range in the parametrization isf")
  }
  ifelse(is.numeric(logPDF), ifelse(log, return(logPDF), return(exp(logPDF))), 
         logPDF)
}


# Cumulative Probability Function

pdtp <- function (x, mu, par1, par2, delta1, delta2, F, f, param = "tp", 
                  log.p = FALSE) 
{
  if (param == "tp") {
    sigma1 = par1
    sigma2 = par2
    eps <- sigma1 * f(0, delta2, log = FALSE)/(sigma1 * 
                                                 f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, 
                                                                                        log = FALSE))
    ifelse(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0, 
           CDF <- ifelse(x < mu, 2 * eps * F((x - mu)/sigma1, 
                                             delta1, log.p = FALSE), eps + (1 - eps) * (2 * 
                                                                                          F((x - mu)/sigma2, delta2, log.p = FALSE) - 
                                                                                          1)), CDF <- "invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma1 = par1 * (1 + par2)
    sigma2 = par1 * (1 - par2)
    eps <- sigma1 * f(0, delta2, log = FALSE)/(sigma1 * 
                                                 f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, 
                                                                                        log = FALSE))
    ifelse(par1 > 0 & abs(par2) < 1 & delta1 > 0 & delta2 > 
             0, CDF <- ifelse(x < mu, 2 * eps * F((x - mu)/sigma1, 
                                                  delta1, log.p = FALSE), eps + (1 - eps) * (2 * F((x - 
                                                                                                      mu)/sigma2, delta2, log.p = FALSE) - 1)), CDF <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma1 = par1 * par2
    sigma2 = par1/par2
    eps <- sigma1 * f(0, delta2, log = FALSE)/(sigma1 * 
                                                 f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, 
                                                                                        log = FALSE))
    ifelse(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0, 
           CDF <- ifelse(x < mu, 2 * eps * F((x - mu)/sigma1, 
                                             delta1, log.p = FALSE), eps + (1 - eps) * (2 * 
                                                                                          F((x - mu)/sigma2, delta2, log.p = FALSE) - 
                                                                                          1)), CDF <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  ifelse(is.numeric(CDF), ifelse(log.p, return(log(CDF)), 
                                 return(CDF)), CDF)
}


# Random Number Generation

rdtp <- function (n, mu, par1, par2, delta1, delta2, rF, f, param = "tp") 
{
  if (param == "tp") {
    sigma1 = par1
    sigma2 = par2
    eps <- sigma1 * f(0, delta2, log = FALSE)/(sigma1 * 
                                                 f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, 
                                                                                        log = FALSE))
    ifelse(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0, 
           sample <- ifelse(runif(n) < eps, mu - sigma1 * abs(rF(n, 
                                                                 delta1)), mu + sigma2 * abs(rF(n, delta2))), 
           sample <- "invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma1 = par1 * (1 + par2)
    sigma2 = par1 * (1 - par2)
    eps <- sigma1 * f(0, delta2, log = FALSE)/(sigma1 * 
                                                 f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, 
                                                                                        log = FALSE))
    ifelse(par1 > 0 & abs(par2) < 1 & delta1 > 0 & delta2 > 
             0, sample <- ifelse(runif(n) < eps, mu - sigma1 * 
                                   abs(rF(n, delta1)), mu + sigma2 * abs(rF(n, delta2))), 
           sample <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma1 = par1 * par2
    sigma2 = par1/par2
    eps <- sigma1 * f(0, delta2, log = FALSE)/(sigma1 * 
                                                 f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, 
                                                                                        log = FALSE))
    ifelse(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0, 
           sample <- ifelse(runif(n) < eps, mu - sigma1 * abs(rF(n, 
                                                                 delta1)), mu + sigma2 * abs(rF(n, delta2))), 
           sample <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  return(sample)
}



# Quantile Function

qdtp <- function (p, mu, par1, par2, delta1, delta2, qF, f, param = "tp") 
{
  if (param == "tp") {
    sigma1 = par1
    sigma2 = par2
    eps <- sigma1 * f(0, delta2, log = FALSE)/(sigma1 * 
                                                 f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, 
                                                                                        log = FALSE))
    ifelse(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0, 
           Q <- ifelse(p < eps, mu + sigma1 * qF(p/(2 * eps), 
                                                 delta1), mu + sigma2 * qF(0.5 * (1 + (p - eps)/(1 - 
                                                                                                   eps)), delta2)), Q <- "invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp")
  }
  if (param == "eps") {
    sigma1 = par1 * (1 + par2)
    sigma2 = par1 * (1 - par2)
    eps <- sigma1 * f(0, delta2, log = FALSE)/(sigma1 * 
                                                 f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, 
                                                                                        log = FALSE))
    ifelse(par1 > 0 & abs(par2) < 1 & delta1 > 0 & delta2 > 
             0, Q <- ifelse(p < eps, mu + sigma1 * qF(p/(2 * 
                                                           eps), delta1), mu + sigma2 * qF(0.5 * (1 + (p - 
                                                                                                         eps)/(1 - eps)), delta2)), Q <- "invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps")
  }
  if (param == "isf") {
    sigma1 = par1 * par2
    sigma2 = par1/par2
    eps <- sigma1 * f(0, delta2, log = FALSE)/(sigma1 * 
                                                 f(0, delta2, log = FALSE) + sigma2 * f(0, delta1, 
                                                                                        log = FALSE))
    ifelse(par1 > 0 & par2 > 0 & delta1 > 0 & delta2 > 0, 
           Q <- ifelse(p < eps, mu + sigma1 * qF(p/(2 * eps), 
                                                 delta1), mu + sigma2 * qF(0.5 * (1 + (p - eps)/(1 - 
                                                                                                   eps)), delta2)), Q <- "invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf")
  }
  return(Q)
}
