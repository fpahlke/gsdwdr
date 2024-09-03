#' 
#' @title Get Double Rejection Design
#' 
#' @description 
#' Function to compute boundaries of double rejection design
#' 
#' @param alpha one-sided type I error rate
#' @param beta type II error rate
#' @param informationRates 2-dimensional vector of information rates
#' @param delayedInformation amount of pipeline information
#' @param typeOfDesign either asOF or asP for alpha-spending
#' @param typeBetaSpending either bsOF or bsP for beta-spending
#' 
#' @return 
#' List containing the boundary set of the input design
#' 
getDoubleRejectionDesign <- function(
        alpha = 0.025,
        beta = 0.2,
        informationRates = c(0.5, 1),
        delayedInformation = 0.3,
        typeOfDesign = c("asOF", "asP"),
        typeBetaSpending = c("bsOF", "bsP")) {
    
    typeOfDesign <- match.arg(typeOfDesign)
    typeBetaSpending <- match.arg(typeBetaSpending)

    # alpha and beta-spending
    aS <- ifelse(typeOfDesign == "asP", pocockSpendingFunction, oBrienFlemingSpendingFunction)
    bS <- ifelse(typeBetaSpending == "bsP", pocockSpendingFunction, oBrienFlemingSpendingFunction)

    # upper boundaries

    # placeholder for upper boundaries
    u <- rep(NA, 2)

    # calculation of d1
    d1 <- qnorm(1 - alpha)

    for (k in 1:2) {
        if (k == 1) {
            informationRatesPlusDelay <- c(informationRates[1], 
                informationRates[1] + delayedInformation)
        } else {
            informationRatesPlusDelay <- informationRates
        }
        u[k] <- uniroot(
            function(x) {
                if (k == 1) {
                    bounds <- matrix(c(
                        x, d1,
                        Inf, Inf
                    ), nrow = 2, byrow = TRUE)
                } else if (k == 2) {
                    bounds <- matrix(c(
                        -Inf, x,
                        u[1], Inf
                    ), nrow = 2, byrow = TRUE)
                }
                probs <- getGroupSequentialProbabilities(bounds, 
                    informationRatesPlusDelay)
                if (k == 1) {
                    # pocockSpendingFunction(Z1 > u1 
                    # AND Z1.tilde > d1) = a(I1), see equation (2) page 7
                    probs[2, 2] - probs[1, 2] - 
                        aS(alpha, informationRates[1])
                } else if (k == 2) {
                    # pocockSpendingFunction(Z1 <= u1 
                    # AND Z2 > d2) = alpha - a(I1), see equation (2) page 7
                    probs[2, 2] - probs[1, 2] - (aS(alpha, informationRates[2]) - 
                            aS(alpha, informationRates[1]))
                }
            },
            lower = -10, upper = 10
        )$root
    }

    # lower boundaries

    # placeholder for lower boundaries
    l1 <- rep(NA, 2)

    # calculation is done via iterative bisection search
    cLower1 <- 0
    cUpper1 <- 100
    precision1 <- 1
    iteration <- 1e5

    while (precision1 > 1e-6) {
        shift <- (cLower1 + cUpper1) / 2
        for (k in 1:2) {
            if (k == 1) {
                informationRatesPlusDelay <- c(informationRates[1], 
                    informationRates[1] + delayedInformation)
            } else {
                informationRatesPlusDelay <- informationRates
            }
            ncp <- matrix(rep(sqrt(informationRatesPlusDelay), 2), 
                nrow = 2, byrow = TRUE) * sqrt(shift)
            precision2 <- 1
            cLower2 <- -8
            cUpper2 <- 8
            while (precision2 > 1e-8) {
                x <- (cLower2 + cUpper2) / 2
                if (k == 1) {
                    bounds <- matrix(c(
                        u[1], d1,
                        Inf, Inf
                    ), nrow = 2, byrow = TRUE) - ncp
                    probs <- getGroupSequentialProbabilities(bounds, 
                        informationRatesPlusDelay)
                    # pnorm(x - ncp[1]) + probs[1, 2] = P_H1(Z1 < l1) + P_H1(Z1 > u1 
                    # AND Z1.tilde <= d1), see equation (2) page 7
                    ifelse(pnorm(x - ncp[1]) + probs[1, 2] < bS(beta, 
                            informationRates[1]),
                        cLower2 <- x, cUpper2 <- x
                    )
                } else if (k == 2) {
                    bounds <- matrix(c(
                        pmin(l1[1], u[1]), x,
                        u[1], Inf
                    ), nrow = 2, byrow = TRUE) - ncp
                    probs <- getGroupSequentialProbabilities(bounds, 
                        informationRates)
                    #  probs[1, 2] = P_H1(Z1 in [l1,u1], Z1.tilde <= d2)
                    ifelse(probs[1, 2] < bS(beta, informationRates[2]) -
                        bS(beta, informationRates[1]),
                    cLower2 <- x, cUpper2 <- x
                    )
                }
                iteration <- iteration - 1
                ifelse(iteration > 0, precision2 <- cUpper2 - cLower2, precision2 <- 0)
            }
            l1[k] <- x
        }
        ifelse(l1[2] < u[2], cLower1 <- shift, cUpper1 <- shift)
        ifelse(iteration > 0, precision1 <- cUpper1 - cLower1, precision1 <- 0)
    }

    design <- list()
    design$alpha <- alpha
    design$beta <- beta
    design$alphaSpending <- alphaSpending
    design$betaSpending <- betaSpending
    design$informationRates <- informationRates
    design$delayedInformation <- delayedInformation
    design$upper <- c(u[1], d1, u[2])
    design$lower <- l1
    design$shift <- shift
    class(design) <- c("DoubleRejectionApproach")
    return(design)
}
