# Function to compute performance characteristics for the double rejection approach
# Inputs:
#   design: design object; output of the function getDoubleRejectionDesign()
#   alternative: vector of effect sizes
#   stDev: standard deviation
#   maxNumberOfSubjects: maximum number of subjects for both groups
# Output:
#   Data frame containing computed performance values for each alternative effect size
getPowerMeansDoubleRejectionApproach <- function(design, alternative, stDev, maxNumberOfSubjects) {
  
  # Number of alternative effect sizes
  nAlternative <- length(alternative)
  
  # Initialize results data frame
  results <- data.frame(
    alternative = alternative,
    rejInterim = rep(NA, nAlternative),
    power = rep(NA, nAlternative),
    futility = rep(NA, nAlternative),
    expectedSampleSize = rep(NA, nAlternative)
  )
  
  # Extract design parameters
  u1 <- design$upper[1]
  d1 <- design$upper[2]
  d2 <- design$upper[3]
  l1 <- design$lower[1]
  
  # Compute covariance matrix from information rates
  informationRatesWithDelay <- c(design$informationRates[1], design$informationRates[1] + design$delayedInformation, 1)
  sigma <- getCovarianceFromInformation(informationRatesWithDelay)
  
  # Function to calculate values for a given effect size
  calculatePerformance <- function(effect) {
    # Compute mean vector
    mu <- effect / sqrt(2 * stDev^2) * sqrt((informationRatesWithDelay * maxNumberOfSubjects / 2))
    
    # Compute interim rejection probability
    rejInterim <- sadmvn(lower = c(u1, d1), 
                    upper = c(Inf, Inf), 
                    mean = mu[1:2], 
                    varcov = sigma[c(1,2), c(1,2)],
                    maxpts = 2000 * d, 
                    abseps = 1e-06, 
                    releps = 0)
    
    # Compute final rejection probability
    rejFinal <- sadmvn(lower = c(l1, d2),
                  upper = c(u1, Inf), 
                  mean = mu[c(1, 3)],
                  varcov = sigma[c(1, 3), c(1, 3)], 
                  maxpts = 2000 * d,
                  abseps = 1e-06, 
                  releps = 0)
    
    # Compute power
    power <- rejInterim + rejFinal
    
    # Compute futility
    futility <- sadmvn(
      lower = c(u1, -Inf),
      upper = c(Inf, d1), 
      mean = mu[1:2], 
      varcov = sigma[c(1,2), c(1,2)],
      maxpts = 2000 * d, 
      abseps = 1e-06,
      releps = 0
    ) + pnorm(l1, mean = mu[1], sd = 1)
    
    # Compute continuation region probability
    continue <- pnorm(u1, mean = mu[1], sd = 1) - pnorm(l1, mean = mu[1], sd = 1)
    
    # Compute expected sample size
    expectedSampleSize <- (rejInterim + futility) * informationRatesWithDelay[2] * maxNumberOfSubjects + continue * maxNumberOfSubjects
    
    return(c(rejInterim, rejFinal, power, futility, expectedSampleSize))
  }
  
  # Apply calculate_values function to each alternative effect size
  results <- t(apply(matrix(alternative, nrow = nAlternative, byrow = TRUE), 1, calculatePerformance))
  
  # Set column names for the result data frame
  colnames(results) <- c("rejInterim", "rejFinal", "power", "futility", "expectedSampleSize")
  
  return(as.data.frame(results))
}
