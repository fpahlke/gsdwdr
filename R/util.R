# O'Brien-Fleming spending function, X is either alpha or beta, I1 <= 1
OBF <- function(X, interimInformation) {
  res <- ifelse(interimInformation == 0, 0, min(2 * (1 - pnorm(qnorm(1 - X / 2) / sqrt(interimInformation))), X))

  return(res)
}

# Pocock spending funtion, X is either alpha or beta, I1 <= 1
P <- function(X, interimInformation) {
  return(min(X * (log(1 + (exp(1) - 1) * interimInformation)), X))
}


# Calculates covariance matrix from information fraction vector
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

# Calculates the probability to obtain futility in the DRGSD design
# Inputs:
#   fut: lower bounds of continuation region, e.g. from a DRGSD object (getDesignGroupSequential(...,bindingFutility = F))
#   crit: upper bounds of continuation region, e.g. from a DRGSD object (getDesignGroupSequential(...,bindingFutility = F))
#   dec: decision critical values, e.g. from a DRGSD object (getDesignGroupSequential(...,bindingFutility = F))
#   alternative: vector of effect values
#   stDev: standard deviation assumption
#   maxNumberOfSubjects: maximum number of subjects for both groups
#   informationRates: 2-dimensional vector of information rates
# Output:
#   Probability to obtain futility given inputs
calcFutility <- function(fut, crit, dec, alternative, stDev, maxNumberOfSubjects, informationRates) {
  probs <- length(alternative)

  sigma <- getCovarianceFromInformation(informationRates)

  for (i in 1:length(alternative)) {
    mu <- alternative[i] / sqrt(2 * stDev^2) * sqrt(informationRates * maxNumberOfSubjects / 2)

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
