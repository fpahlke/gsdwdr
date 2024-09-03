#'
#' O'Brien-Fleming spending function
#' 
#' @param x is either alpha or beta, I1 <= 1
#' 
oBrienFlemingSpendingFunction <- function(x, interimInformation) {
  res <- ifelse(interimInformation == 0, 
      0, 
      min(2 * (1 - pnorm(qnorm(1 - x / 2) / sqrt(interimInformation))), x)
  )
  return(res)
}

#'
#' Pocock spending function
#' 
#' @param x is either alpha or beta, I1 <= 1
#' 
pocockSpendingFunction <- function(x, interimInformation) {
  return(min(x * (log(1 + (exp(1) - 1) * interimInformation)), x))
}


#' 
#' Calculates covariance matrix from information fraction vector
#' 
getCovarianceFromInformation <- function(information) {
  sigma <- matrix(NA, ncol = length(information), nrow = length(information))

  for (k in seq_along(information)) {
    for (h in seq_along(information)) {
      sigma[k, h] <- sqrt(information[k] / information[h])

      if (sigma[k, h] > 1) sigma[k, h] <- 1 / sigma[k, h]
    }
  }

  return(sigma)
}

#' 
#' @title Calc Futility
#' 
#' @description 
#' Calculates the probability to obtain futility in the DRGSD design
#' 
#' @param fut: lower bounds of continuation region, e.g. from a DRGSD object (getDesignGroupSequential(..., bindingFutility = FALSE))
#' @param crit: upper bounds of continuation region, e.g. from a DRGSD object (getDesignGroupSequential(..., bindingFutility = FALSE))
#' @param dec: decision critical values, e.g. from a DRGSD object (getDesignGroupSequential(..., bindingFutility = FALSE))
#' @param alternative: vector of effect values
#' @param stDev: standard deviation assumption
#' @param maxNumberOfSubjects: maximum number of subjects for both groups
#' @param informationRates: 2-dimensional vector of information rates
#' 
#' @return 
#' Probability to obtain futility given inputs
#' 
calcFutility <- function(fut, crit, dec, alternative, stDev, 
    maxNumberOfSubjects, informationRates) {
  probs <- length(alternative)

  sigma <- getCovarianceFromInformation(informationRates)

  for (i in 1:length(alternative)) {
    mu <- alternative[i] / sqrt(2 * stDev^2) * 
        sqrt(informationRates * maxNumberOfSubjects / 2)

    probs[i] <- sadmvn(
      lower = c(-Inf, -Inf),
      upper = c(fut[1], dec[1]),
      mean = mu[c(1, 2)],
      varcov = sigma[c(1, 2), c(1, 2)],
      maxpts = 2000 * d, abseps = 1e-06, releps = 0
    ) +
      sadmvn(
        lower = c(crit[1], -Inf),
        upper = c(Inf, dec[1]),
        mean = mu[c(1, 2)],
        varcov = sigma[c(1, 2), c(1, 2)],
        maxpts = 2000 * d, abseps = 1e-06, releps = 0
      )
  }

  return(probs)
}
