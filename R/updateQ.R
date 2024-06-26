#' Update transition rates
#'
#' @param nbStates Number of states
#' @param data Matrix of all data (observations and transitions)
#' @param switch Matrix of transitions
#' @param priorShape Shape of prior gamma distribution for the transition rates
#' @param priorRate Rate of prior gamma distribution for the transition rates
#' @param priorCon Concentration of prior Dirichlet distribution for
#' transition probabilities
#'
#' @return Updated transition matrix
#'
#' @importFrom gtools rdirichlet
#' @export
updateQ <- function(nbStates, data, switch, priorShape, priorRate, priorCon, kappa) {
  # time spent in each state
  stateTab <- rbind(
    data[1, c("time", "state")], # to include first interval
    switch,
    data[nrow(data), c("time", "state")]
  ) # to include last interval
  timeInStates <- sapply(1:nbStates, function(i) {
    sum(diff(stateTab[, "time"])[which(stateTab[, "state"] == i)], na.rm = TRUE)
  })

  # count intervals spent in each state
  countIntervals <- table(factor(rle(data[, "state"])$values, levels = 1:nbStates))

  # sample rates out of each state
  shape <- priorShape + countIntervals
  rate <- priorRate + timeInStates
  r <- pmin(rgamma(n = nbStates, shape = shape, rate = rate), kappa)

  # sample transition probabilities
  allCounts <- table(factor(data[-nrow(data), "state"], levels = 1:nbStates), factor(data[-1, "state"], levels = 1:nbStates))
  nonDiagCounts <- matrix(t(allCounts)[!diag(nbStates)], nrow = nbStates, byrow = TRUE)
  trProbs <- t(sapply(seq_len(nrow(nonDiagCounts)), function(i) rdirichlet(n = 1, alpha = nonDiagCounts[i, ] + priorCon[-i])))
  # there is an NaN issue when any states are not observed
  # we can either set these all equal to 1/(nbStates-1) or I think more realistically,
  # randomly assign all the probability to one transition
  for (i in which(is.nan(rowSums(trProbs)))) {
    trProbs[i, ] <- sample(c(1, rep(0, nbStates - 2)), nbStates - 1, FALSE)
  }

  # update generator matrix from rates and transition probs
  Q <- -diag(nbStates)
  Q[!Q] <- t(trProbs)
  Q <- t(Q * rep(r, each = nbStates))

  return(Q)
}
