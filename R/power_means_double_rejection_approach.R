#'
#' @title Get Power Means Double Rejection Approach
#'
#' @description
#' Function to compute performance characteristics for the double rejection approach
#'
#' @param design design object; output of the function getDoubleRejectionDesign()
#' @param alternative vector of effect sizes
#' @param stDev standard deviation
#' @param maxNumberOfSubjects maximum number of subjects for both groups
#'
#' @return
#' Data frame containing computed performance values for each alternative effect size
#'
getPowerMeansDoubleRejectionApproach <- function(
        design, alternative, stDev, maxNumberOfSubjects) {
        
    # number of alternative effect sizes
    nAlternative <- length(alternative)

    # initialize results data frame
    results <- data.frame(
        alternative = alternative,
        rejInterim = rep(NA, nAlternative),
        power = rep(NA, nAlternative),
        futility = rep(NA, nAlternative),
        expectedSampleSize = rep(NA, nAlternative)
    )

    # extract design parameters
    u1 <- design$upper[1]
    d1 <- design$upper[2]
    d2 <- design$upper[3]
    l1 <- design$lower[1]

    # compute covariance matrix from information rates
    informationRatesWithDelay <- c(
        design$informationRates[1],
        design$informationRates[1] + design$delayedInformation, 1
    )
    sigma <- getCovarianceFromInformation(informationRatesWithDelay)

    # function to calculate values for a given effect size
    calculatePerformance <- function(effect) {
        # compute mean vector
        mu <- effect / sqrt(2 * stDev^2) *
            sqrt(informationRatesWithDelay * maxNumberOfSubjects / 2)

        # compute interim rejection probability
        rejInterim <- sadmvn(
            lower = c(u1, d1),
            upper = c(Inf, Inf),
            mean = mu[1:2],
            varcov = sigma[c(1, 2), c(1, 2)],
            maxpts = 2000 * d,
            abseps = 1e-06,
            releps = 0
        )

        # compute final rejection probability
        rejFinal <- sadmvn(
            lower = c(l1, d2),
            upper = c(u1, Inf),
            mean = mu[c(1, 3)],
            varcov = sigma[c(1, 3), c(1, 3)],
            maxpts = 2000 * d,
            abseps = 1e-06,
            releps = 0
        )

        # compute power
        power <- rejInterim + rejFinal

        # compute futility
        futility <- sadmvn(
            lower = c(u1, -Inf),
            upper = c(Inf, d1),
            mean = mu[1:2],
            varcov = sigma[c(1, 2), c(1, 2)],
            maxpts = 2000 * d,
            abseps = 1e-06,
            releps = 0
        ) + pnorm(l1, mean = mu[1], sd = 1)

        # compute continuation region probability
        continue <- pnorm(u1, mean = mu[1], sd = 1) - pnorm(l1, mean = mu[1], sd = 1)

        # compute expected sample size
        expectedSampleSize <- (rejInterim + futility) * informationRatesWithDelay[2] *
            maxNumberOfSubjects + continue * maxNumberOfSubjects

        return(c(rejInterim, rejFinal, power, futility, expectedSampleSize))
    }

    # apply calculate_values function to each alternative effect size
    results <- t(apply(matrix(alternative,
        nrow = nAlternative,
        byrow = TRUE
    ), 1, calculatePerformance))

    # set column names for the result data frame
    colnames(results) <- c(
        "rejInterim",
        "rejFinal",
        "power",
        "futility",
        "expectedSampleSize"
    )

    return(as.data.frame(results))
}
