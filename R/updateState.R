#' Update state sequence
#'
#' @param obs Matrix of observations, with columns `x`, `y`, `time`, and `state`
#' @param knownStates Vector of known, fixed states
#' @param switch Matrix of state switches, with columns `time` and `state`
#' @param updateLim Vector of two elements: min and max length of updated interval
#' (each an integer, corresponding to a number of observations)
#' @param updateProbs Vector of probabilities of picking each number from
#' updateLim[1] to updateLim[2], for the length of the update interval
#' @param kappa Upper bounds of transition rate
#'
#' @return List of two elements:
#'   * newSwitch: Updated matrix of transitions
#'   * newData: Updated matrix of all data (observations and transitions)
#'
#' @details The update is done with the function \code{sample_path} from the
#' package ECctmc (Fintzi, 2017).
#'
#' @references
#' Jon Fintzi (2017). ECctmc: Simulation from Endpoint-Conditioned Continuous
#' Time Markov Chains. R package version 0.2.4.
#' https://CRAN.R-project.org/package=ECctmc
#'
#' TODO: needs priors
#' @export
updateState <- function(obs, nbStates, knownStates, switch, updateLim, param, mu, Hmat, updateProbs = NULL, rateparam = NULL, kappa = NULL, model = NULL) {
  nbObs <- nrow(obs)

  if (is.null(updateProbs)) {
    updateProbs <- rep(1, length(updateLim[1]:updateLim[2])) / length(updateLim[1]:updateLim[2])
  }

  # select interval to update
  if (updateLim[1] < updateLim[2]) {
    len <- sample(updateLim[1]:updateLim[2], size = 1, prob = updateProbs)
  } else {
    len <- updateLim[1]
  }
  begin <- sample(1:(nbObs - len), size = 1)
  end <- begin + len
  Tbeg <- obs[begin, "time"]
  Tend <- obs[end, "time"]

  # sample state sequence conditional on start and end state
  path <- sample_path_mr_(
    a = obs[begin, "state"],
    b = obs[end, "state"],
    t0 = Tbeg,
    t1 = Tend,
    lng0 = obs[begin, "x"],
    lat0 = obs[begin, "y"],
    lng1 = obs[end, "x"],
    lat1 = obs[end, "y"],
    group = obs[begin, "group"],
    k = kappa,
    nbStates = nbStates,
    param = param,
    mu = mu,
    Hmat = Hmat[c(begin, end), ],
    rateparam = rateparam,
    model = ifelse(is.na(model), "NA", model)
  )

  # check for autoreject
  if (nrow(path) == 1 && all(path == -1)) {
    stop("Exceeded path simulation iteration limit")
  }

  path <- path[-c(1, nrow(path)), ] # remove 1st and last rows (observations)

  # update state sequence
  newSwitch <- as.matrix(
    rbind(
      switch[switch[, "time"
      ] < Tbeg,
      ],
      path,
      switch[switch[, "time"
      ] > Tend,
      ],
      deparse.level = 0
    ),
    ncol = 2, drop = FALSE
  )

  # remove switches into same state
  fakeSwitch <- which(newSwitch[-1, "state", drop = FALSE] == newSwitch[-nrow(newSwitch), "state", drop = FALSE]) + 1

  if (length(fakeSwitch) > 0) {
    newSwitch <- newSwitch[-fakeSwitch, , drop = FALSE]
  }
  if (nrow(newSwitch) > 0) {
    newData <- rbind(
      obs,
      cbind("x" = NA, "y" = NA, "time" = newSwitch[, "time"], "ID" = rep(obs[1, "ID"], nrow(newSwitch)), "state" = newSwitch[, "state"], "group" = rep(obs[1, "group"], nrow(newSwitch)))
    )
    rownames(newData) <- NULL
  } else {
    newData <- obs
  }

  newData <- newData[order(newData[, "time"]), c("x", "y", "time", "ID", "state", "group")]

  # update state sequence for new switches
  ind <- which(newData[, "time"] > Tbeg & newData[, "time"] < Tend)
  for (t in ind) {
    if (!is.na(newData[t, "x"])) {
      newData[t, "state"] <- newData[t - 1, "state"]
    }
  }

  # knownStates override
  newData[which(!is.na(newData[, "x"])), ][which(!is.na(knownStates)), "state"] <- knownStates[which(!is.na(knownStates))]

  return(list(newSwitch = newSwitch, newData = newData, len = len))
}
