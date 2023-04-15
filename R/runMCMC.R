
#' Run MCMC iterations
#'
#' @param track data.frame. Data, with columns \code{"x"}, \code{"y"}, \code{"time"}, and \code{"ID"}
#' @param nbStates integer. Number of states
#' @param nbIter integer. Number of iterations for the MCMC
#' @param inits list. Initial parameters:
#' \itemize{
#'   \item{"tau_pos":} vector. Initial tau_pos for each state, of length \code{"nbStates"}
#'   \item{"tau_vel":} vector. Initial tau_vel for each state, of length \code{"nbStates"}
#'   \item{"sigma":} vector. Initial sigma for each state, of length \code{"nbStates"}
#'   \item{"mu":} list. Initial OUF range centre coordinate pairs \code{"(x, y)"}, with \code{"NA"} for IOU states
#'   \item{"Q":}
#'   \item{"state":} vector. Initial state sequence, length \code{"nrow(track)"}
#'   \item{"alpha":}
#'   \item{"t_alpha":}
#' }
#' @param fixed list of fixed parameters:
#' \itemize{
#'   \item{"tau_pos":} vector. Fixed values for tau_pos for each state with \code{"NA"} for parameters to estimate. Length \code{"nbStates"}
#'   \item{"tau_vel":} vector. Fixed values for tau_vel for each state with \code{"NA"} for parameters to estimate. Length \code{"nbStates"}
#'   \item{"sigma":} vector. Fixed values for sigma for each state with \code{"NA"} for parameters to estimate. Length \code{"nbStates"}
#'   \item{"mu":} list. Fixed OUF range centre coordinate pairs \code{"(x, y)"}, with \code{"NA"} for IOU states or pairs to estimate
#'   \item{"Q":} Unused
#'   \item{"knownStates":} vector. Known states with \code{"NA"} for those to estimate. Length \code{"nrow(track)"}
#'   \item{"kappa":}
#' }
#' @param priors list. Parameters of prior distributions, with components:
#' \itemize{
#'   \item{"mean":} vector. Means for normal priors on movement parameters, of length \code{"5*nbStates"}
#'   \item{"sd":} vector. Standard deviations for normal priors on movement parameters, of length \code{"5*nbStates"}
#'   \item{"shape":} vector. Shapes of gamma priors for the transition rates
#'   \item{"rate":} vector. Rates of gamma priors for the transition rates
#'   \item{"con":} vector. Concentrations of Dirichlet priors for transition probabilities
#'   \item{"rate_mean":}
#'   \item{"rate_sd":}
#' }
#' @param props list. Parameters of proposal distributions, with components:
#' \itemize{
#'   \item{"S":} matrix. Initial value for the lower triangular matrix of RAM algorithm, so that the covariance matrix of the proposal distribution is \code{"SS'"}.
#'   \item{"updateLim":} vector. Two values: min and max length of updated state sequence
#'   \item{"updateProbs":} vector. Probabilities for each element of \code{"updateLim[1]:updateLim[2]"} (if \code{"NULL"},
#'   all values are equiprobable)
#'   \item{"rate_S":}
#' }
#' @param tunes list. Tuning parameters, with components:
#' \itemize{
#'   \item{"thinStates":} integer. Thinning factor for the posterior state sequences (needed because
#' of memory limitations)
#' }
#' @param Hmat matrix. Observation error variance (four columns, and one row
#' for each row of data)
#' @param updateState logical. If \code{"FALSE"}, the state process is not updated
#' (for exploratory analysis only, useful for testing single state models)
#' @param adapt integer. If \code{"adapt"} > 0, use the the Robust Adaptive Metropolis (RAM) algorithm
#' by Vihola (2012) to update the proposal distribution for each parameter at each iteration (up to \code{"adapt"} iterations)
#' to a target acceptance rate of 23.4\%.
#' @param model
#'
#' @references
#' Michelot, T., Blackwell, P.G. (2019).
#' State‐switching continuous‐time correlated random walks.
#' Methods Ecol Evol, 10: 637-649. doi:10.1111/2041-210X.13154
#'
#' Vihola, M. (2012).
#' Robust adaptive Metropolis algorithm with coerced acceptance rate.
#' Stat Comput, 22: 997-1008. doi:10.1007/s11222-011-9269-5
#' @examples
#' \dontrun{
#' }
#'
#' @importFrom stats dnorm runif rnorm rexp rgamma
#' @importFrom prettyunits vague_dt pretty_dt
#' @importFrom ramcmc adapt_S
#' @export
#'
#' @useDynLib MSctmm
runMCMC <- function(track, nbStates, nbIter, inits, fixed, priors,
                    props, tunes, Hmat,
                    updateState = TRUE, adapt = FALSE, model = NA)
{
  # Check track df
  if (!is.data.frame(track))
    stop("argument 'track' is not a data.frame")
  if (any(!(c("x", "y", "time") %in% colnames(track))))
    stop("argument 'track' is missing required column(s): ", paste(c("x", "y", "time")[which(!(c("x", "y", "time") %in% colnames(track)))], collapse = ", "))

  if (!("ID" %in% colnames(track))) {
    warning("argument 'track' should have column 'ID', assuming data is one individual")
    track$ID <- 1
  } else {
    track$ID <- as.numeric(as.factor(track$ID))
  }

  # Check inits arguments and lengths
  if (is.null(inits$tau_pos) | is.null(inits$tau_vel) | is.null(inits$sigma) | is.null(inits$mu) | is.null(inits$state))
    stop("argument 'inits' missing: ", paste(c("tau_pos", "tau_vel", "sigma", "mu", "state")[which(sapply(inits[c("tau_pos", "tau_vel", "sigma", "mu", "state")], is.null))], collapse = ", "))
  for (arg in c("tau_pos", "tau_vel", "sigma")) {
    if (length(inits[[arg]]) != nbStates)
      stop("argument 'inits$", arg,"' has the wrong length, expected ", nbStates, ", but got ", length(inits[[arg]]))
  }
  if (all(sapply(inits$mu, length) != rep(nbStates, 2))) {
    stop("argument 'inits$mu' has the wrong length, expected ", paste(rep(nbStates, 2), collapse = " "), ", but got", paste(sapply(inits$mu, length), collapse = " "))
  }
  if (length(inits$state) != nrow(track))
    stop("'inits$state' has the wrong length, expected ", nrow(track), ", but got ", length(inits$state))
  if (is.null(inits$Q) & (is.null(inits$kappa) | is.null(inits$alpha) | is.null(inits$t_aplha))) {
    stop("argument 'inits$Q' is null, expected ", paste(paste0("inits$", c("kappa", "alpha", "t_alpha")[which(sapply(inits[c("kappa", "alpha", "t_alpha")], is.null))]), collapse = ", "), " to be specified")
  }

  # TODO: Check Q length/dims in both init and fixed
  # TODO: Check not Q AND model NULL/NA one or the other
  # TODO: Use missing(model) instead?

  # Check fixed arguments and lengths
  for (arg in c("tau_pos", "tau_vel", "sigma")) {
    if (is.null(fixed[[arg]])) {
      fixed[[arg]] <- rep(NA, nbStates)
    }
    if (length(fixed[[arg]]) != nbStates)
      stop("argument 'fixed$", arg,"' has the wrong length, expected ", nbStates, ", but got ", length(fixed[[arg]]))
  }
  if (is.null(fixed$mu)) {
    fixed$mu <- rep(list(c(NA,NA)), nbStates)
  } else {
    if (all(sapply(fixed$mu, length) != rep(nbStates, 2))) {
      stop("argument 'fixed$mu' has the wrong length, expected ", paste(rep(nbStates, 2), collapse = " "), ", but got", paste(sapply(fixed$mu, length), collapse = " "))
    }
  }
  if (is.null(fixed$knownStates)) {
    fixed$knownStates <- rep(NA, nrow(track))
  } else {
    if (length(fixed$knownStates) != nrow(track)) {
      stop("argument 'fixed$knownStates' has the wrong length, expected ", nrow(track), ", but got", length(fixed$knownStates))
    }
  }

  # TODO update this based on length of mean and sd? if 5 * nbStates then shape rate, con needed
  # If mean and sd includes rate_mean rate_sd, and model not null, etc, then shape rate con not needed

  # Check priors arguments and lengths
  if (is.null(priors$mean) | is.null(priors$sd) | is.null(priors$shape) | is.null(priors$rate) | is.null(priors$con))
    stop("argument 'priors' missing: ", paste(c("mean", "sd", "shape", "rate", "con")[which(sapply(priors[c("mean", "sd", "shape", "rate", "con")], is.null))], collapse = ", "))
  for (arg in c("mean", "sd")) {
    if (length(priors[[arg]]) != 5 * nbStates)
      stop("argument 'priors$", arg,"' has the wrong length, expected ", 5 * nbStates, ", but got ", length(priors[[arg]]))
  }
  for (arg in c("shape", "rate")) {
    if (length(priors[[arg]]) != nbStates)
      stop("argument 'priors$", arg,"' has the wrong length, expected ", nbStates, ", but got ", length(priors[[arg]]))
  }
  # TODO: Check con length? 1 or nbStates?

  # TODO: CHECK S dims?

  # Check props arguments and lengths
  if (is.null(props$S) | is.null(props$updateLim))
    stop("argument 'props' missing: ", paste(c("S", "updateLim")[which(sapply(props[c("S", "updateLim")], is.null))], collapse = ", "))
  if ("list" %in% class(props$updateLim) & all(sapply(props$updateLim, length) != rep(2, length(unique(track$ID)))))
    stop("argument 'props$updateLim' has the wrong length, expected ", paste(rep(2, length(unique(track$ID))), collapse = " "), " but got ", paste(sapply(props$updateLim, length), collapse = " "))
  if (!"list" %in% class(props$updateLim) & length(props$updateLim) != 2)
    stop("argument 'props$updateLim' has the wrong length, expected ", 2, " but got ", length(props$updateLim))
  if (!is.null(props$updateProbs) & class(props$updateLim) != class(props$updateProbs))
    stop("argument 'props$updateProbs' provided but of wrong type, expected ", class(props$updateLim), " but got ", class(props$updateProbs))
  if (!is.null(props$updateProbs) & "list" %in% class(props$updateProbs) & all(sapply(props$updateLim, length) != sapply(props$updateProbs, length)))
    stop("argument 'props$updateProbs' has the wrong length, expected ", paste(sapply(props$updateLim, length), collapse = " "), " but got ", paste(sapply(props$updateProbs, length), collapse = " "))
  if (!is.null(props$updateProbs) & !"list" %in% class(props$updateProbs) & length(props$updateProbs) != 2)
    stop("argument 'props$updateProbs' has the wrong length, expected ", 2, " but got ", length(props$updateProbs))

  # TODO warn when ignoring other argument (model present, etc.)

  ######################
  ## Unpack arguments ##
  ######################

  # unpack fixed parameters
  fixPar <- fixed[c("tau_pos", "tau_vel", "sigma")]
  fixMu <- fixed$mu
  knownStates <- fixed$knownStates

  # initial movement parameters
  # NB considering only isotropic case for now, this will eventually need to be a matrix
  # and kalman_rcpp and MakeQ/Sigma will need to accept a matrix
  param <- unlist(inits[c("tau_pos", "tau_vel", "sigma")])
  mu <- unlist(inits$mu)
  state <- inits$state

  # overwrite fixed params
  param[which(!is.na(unlist(fixPar)))] <- unlist(fixPar)[which(!is.na(unlist(fixPar)))]
  mu[which(!is.na(unlist(fixMu)))] <- unlist(fixMu)[which(!is.na(unlist(fixMu)))]
  state[which(!is.na(knownStates))] <- knownStates[which(!is.na(knownStates))]

  # unpack prior parameters
  priorMean <- priors$mean
  priorSD <- priors$sd
  priorShape <- priors$shape
  priorRate <- priors$rate
  priorCon <- priors$con

  # unpack proposal parameters
  S <- props$S
  # TODO: CLEAN THIS UP
  if ("list" %in% class(props$updateLim) & length(props$updateLim) == length(unique(track$ID))) {
    updateLim <- lapply(1:length(unique(track$ID)), function(i) {
      lim <- ceiling(props$updateLim[[i]] * nrow(track[which(track$ID == unique(track$ID)[i]), ]))
      if (lim[1] <= 3) lim = lim + 1
      if (lim[1] == lim[2]) lim[2] = lim[2] + 1
      lim
    })
  } else {
    updateLim <- lapply(unique(track$ID), function(id) {
      lim <- ceiling(props$updateLim * nrow(track[which(track$ID == id), ]))
      if (lim[1] <= 3) lim = lim + 1
      if (lim[1] == lim[2]) lim[2] = lim[2] + 1
      lim
    } )
  }
  if (!is.null(props$updateProbs) & "list" %in% class(props$updateProbs) & length(props$updateProbs) == length(unique(track$ID))) {
    updateProbs <- props$updateProbs
    #} else {
    #  updateProbs <- rep(list(props$updateProbs), length(unique(track$ID)))
  } else if (is.null(props$updateProbs)) {
    updateProbs <- lapply(1:length(unique(track$ID)), function(i) {
      rep(1,length(updateLim[[i]][1]:updateLim[[i]][2]))/length(updateLim[[i]][1]:updateLim[[i]][2])
    })
  }

  # unpack tuning parameters
  thinStates <- ifelse(is.null(tunes$thinStates), 1, tunes$thinStates)

  ####################
  ## Initialisation ##
  ####################

  # Prepare data structures
  ids <- unique(track$ID)
  obs <- list()
  known <- list()
  switch <- list()
  data.list <- list()

  # TODO incorporate rate props and priors into main vectors and check them with above
  # check if rate matrix provided
  if (!is.null(inits$Q) & length(inits$Q) == length(ids) & "list" %in% class(inits$Q)) {
    Q <- inits$Q
    names(Q) <- ids
    kappa <- NULL
    rate_priorMean <- NULL
    rate_priorSD <- NULL
    rateparam <- NULL
  } else if (!is.null(inits$Q) & is.null(fixed$kappa)) {
    Q <- rep( list( inits$Q ), length( unique( track$ID ) ) )
    names(Q) <- ids
    kappa <- NULL
    rate_priorMean <- NULL
    rate_priorSD <- NULL
    rateparam <- NULL
  } else { # if no rate matrix provided, rate params must be
    Q <- NULL
    kappa <- fixed$kappa
    rate_S <- props$rate_S
    rate_priorMean <- priors$rate_mean
    rate_priorSD <- priors$rate_sd
    rateparam <- c(inits$alpha, inits$t_alpha)
  }

  names(updateLim) <- ids
  if (exists("updateProbs")) names(updateProbs) <- ids

  for (id in ids) {
    nbObs <- nrow(track[which(track$ID == id), ])
    obs[[id]] <- cbind(
      "x" = track[which(track$ID == id), "x"],
      "y" = track[which(track$ID == id), "y"],
      "time" = track[which(track$ID == id), "time"],
      "ID" = as.numeric(track[which(track$ID == id), "ID"]),
      "state" = state[which(track$ID == id)])
    known[[id]] <- knownStates[which(track$ID == id)]
    #colnames(obs[[id]]) <- c("x", "y", "time", "ID",  "state")
    indSwitch <- which(obs[[id]][-1, "state"] != obs[[id]][-nbObs, "state"]) + 1
    switch[[id]] <- cbind(
      "time" = obs[[id]][indSwitch, "time"] - 0.001,
      "state" = rle(obs[[id]][ , "state"])$values[-1])
    #colnames(switch[[id]]) <- c("time", "state")

    if (!all(is.na(switch[[id]]))) {
      data.list[[id]] <- rbind(
        obs[[id]][ , c( "x", "y", "time", "ID",  "state")],
        cbind(
          "x" = NA,
          "y" = NA,
          "time" = switch[[id]][ , "time"],
          "ID" = id,
          "state" = switch[[id]][ , "state"]
        )
      )
      data.list[[id]] <- data.list[[id]][order(data.list[[id]][ , "time"]), ]
    } else {
      data.list[[id]] <- obs[[id]][ , c("x", "y", "time", "ID",  "state")]
    }
  }

  # flatten data
  data.mat <- do.call("rbind", data.list)
  data.mat <- data.mat[ , c("x", "y", "time", "ID", "state")]

  # initialise Hmat (rows of 0s for transitions)
  HmatAll <- matrix(0, nrow(data.mat), 4)
  HmatAll[which(!is.na(data.mat[ , "x"])), ] <- Hmat

  # initial likelihood
  oldllk <- kalman_rcpp(data = data.mat, param = param, fixmu = mu, Hmat = HmatAll)$llk

  # initial log-prior
  oldlogprior <- getLogPrior(
    param, mu, fixPar, fixMu, priorMean, priorSD,
    rateparam, rate_priorMean, rate_priorSD, kappa, model
  )

  ###############################
  ## Loop over MCMC iterations ##
  ###############################

  allParam <- matrix(NA, nrow = nbIter, ncol = 5 * nbStates)
  colnames(allParam) <- c(
    paste("tau_pos[", 1:nbStates, "]", sep = ""),
    paste("tau_vel[", 1:nbStates, "]", sep = ""),
    paste("sigma[", 1:nbStates, "]", sep = ""),
    paste(c("mu_x[", "mu_y["), rep(1:nbStates, each = 2), c("]", "]"), sep = "")
  )
  allLLk <- rep(NA, nbIter)

  accParam <- rep(0, nbIter)
  if (updateState) {
    allStates <- matrix(NA, nrow = nbIter / thinStates, ncol = nrow(track)) # uses a lot of memory if not thinning!
    allLen <- matrix(NA, nrow = nbIter, ncol = length(ids))
    accSwitch <- rep(0, nbIter)
    if (!is.null(Q)) {
      allRates <- array(NA, dim = c(nbIter, nbStates * (nbStates - 1), length(ids)))
    }
  }


  timing <- matrix(NA, nrow = nbIter / thinStates, ncol = 2)

  t0 <- Sys.time()
  for (iter in 1:nbIter) {
    if (iter %% 100 == 0) {
      cat(
        "\rIteration ", iter, "/", nbIter, "... ",
        vague_dt(difftime(Sys.time(), t0, units = "secs") / iter * (nbIter - iter), "short"),
        " remaining (est)",
        ifelse(
          accSwitch,
          paste0(" -- accSwitch = ", round(sum(accSwitch) / iter * 100), "%"),
          ""
        ),
        " -- accParam = ", round(sum(accParam) / iter * 100), "%",
        "          ", sep = ""
      )
    }


    ################################################
    ## 1. Update discrete state process and rates ##
    ################################################

    if (updateState) {

      newData.list <- data.list
      newSwitch <- switch

      if (is.null(Q) & !is.na(model)) {
        # propose new rate params
        # On working scale [-Inf, Inf]
        rate_u <- rnorm(length(rateparam))
        rate_thetas <- log(rateparam) + as.vector(rate_S %*% rate_u)
        # On natural scale [0, Inf]
        rate_thetasprime <- exp(rate_thetas)
      } else {
        rate_thetas <- NULL
        rate_thetasprime <- rate_thetas
      }

      for (id in ids) {
        upState <- updateState(
          obs = obs[[id]],
          nbStates = nbStates,
          knownStates = known[[id]],
          switch = switch[[id]],
          updateLim = updateLim[[id]],
          updateProbs = updateProbs[[id]],
          Q = ifelse(is.null(Q), NULL, Q[[id]]),
          rateparam = rate_thetasprime,
          kappa = kappa,
          model = model
        )

        newData.list[[id]] <- upState$newData
        newSwitch[[id]] <- upState$newSwitch
        allLen[iter, which(ids == id)] <- upState$len
      }

      # flatten data
      newData.mat <- do.call("rbind", newData.list)

      # update Hmat (rows of 0s for transitions)
      newHmatAll <- matrix(0, nrow(newData.mat), 4)
      newHmatAll[which(!is.na(newData.mat[ , "x"])), ] <- Hmat

      # Calculate acceptance ratio
      newllk <- kalman_rcpp(data = newData.mat, param = param, fixmu = unlist(mu), Hmat = newHmatAll)$llk
      newlogprior <- getLogPrior(
        param, mu, fixPar, fixMu, priorMean, priorSD,
        rate_thetas, rate_priorMean, rate_priorSD, kappa, model
      )
      logHR <- newllk + newlogprior - oldllk - oldlogprior

      if (log(runif(1)) < logHR) {
        # Accept new state sequence
        accSwitch[iter] <- 1
        switch <- newSwitch
        data.list <- newData.list
        obs <- lapply(data.list, function(data) { data[!is.na(data[ , "x"]), ] })
        oldllk <- newllk
        oldlogprior <- newlogprior
        HmatAll <- newHmatAll
      }

      if (!is.null(Q) & is.na(model)) {
        Q <- lapply(ids, function(id) {
          updateQ(
            nbStates = nbStates, data = data.list[[id]], switch = switch[[id]],
            priorShape = priorShape, priorRate = priorRate,
            priorCon = priorCon
          )
        })
        names(Q) <- ids
      }

      #if (adapt & iter >= 1000 & iter <= adapt) {
      if (adapt & iter > 1 & iter <= adapt) {
        rate_S <- adapt_S(rate_S, rate_u, min(1, exp(logHR)), iter)
      }
    }

    ###################################
    ## 2. Update movement parameters ##
    ###################################

    # TODO update mu and include in logprior? make sure correctly bounded?

    pass = F
    while (!pass) {  # ensure tau_p >= tau_v
      # On working scale [-Inf,Inf]
      param_u <- rnorm(length(param))
      param_thetas <- log(param) + as.vector(S %*% param_u)
      # On natural scale [0, Inf]
      param_thetasprime <- unlist(fixPar)
      param_thetasprime[is.na(unlist(fixPar))] <- exp(param_thetas[is.na(unlist(fixPar))])

      # is there potential to get stuck here?
      if (all(param_thetasprime[1:nbStates] >= param_thetasprime[(nbStates + 1):(2 * nbStates)])) {
        pass <- T
      }
    }

    # Calculate acceptance ratio
    data.mat <- do.call("rbind", data.list)
    data.mat <- data.mat[ , c("x", "y", "time", "ID", "state")]
    newlogprior <- getLogPrior(
      param_thetas, mu, fixPar, fixMu, priorMean, priorSD,
      rateparam, rate_priorMean, rate_priorSD, kappa, model
    )
    kalman <- kalman_rcpp(data = data.mat, param = param_thetasprime, fixmu = mu, Hmat = HmatAll)
    newllk <- kalman$llk
    mu <- as.vector(t(kalman$mu))
    logHR <- newllk + newlogprior - oldllk - oldlogprior

    if (log(runif(1)) < logHR) {
      # Accept new parameter values
      accParam[iter] <- 1
      param <- param_thetasprime
      oldllk <- newllk
      oldlogprior <- newlogprior
    }

    #if (adapt & iter >= 1000 & iter <= adapt) {
    if (adapt & iter <= adapt) {
      #S[is.na(unlist(fixPar)), is.na(unlist(fixPar))] <- adapt_S(S[is.na(unlist(fixPar)), is.na(unlist(fixPar))], param_u[is.na(unlist(fixPar))], min(1, exp(logHR)), iter)
      # calculate S by state instead
      for (i in 1:nbStates) {
        index <- seq(i, length(param), by = nbStates)
        S[index,index][is.na(unlist(fixPar)[index]), is.na(unlist(fixPar)[index])] <- adapt_S(S[index,index][is.na(unlist(fixPar)[index]), is.na(unlist(fixPar)[index])], param_u[index][is.na(unlist(fixPar)[index])], min(1, exp(logHR)), iter)
      }
    }

    #########################
    ## Save posterior draw ##
    #########################
    allParam[iter,] <- cbind(matrix(param, ncol = 3 * nbStates), matrix(mu, ncol =  2 * nbStates))
    if (!is.null(Q)) {
      allRates[iter, , ] <- matrix(unlist(lapply(Q, function(q){ q[!diag(nbStates)]})), ncol = length(ids), nrow = nbStates * (nbStates - 1))
    } else {
      allRateParam[iter, ] <- rateparam
    }

    if (iter %% thinStates == 0) {
      allStates[iter / thinStates, ] <- unlist(lapply(obs, function(ob) { ob[ , "state"] }))
      timing[iter / thinStates, ] <- c(iter, Sys.time())
    }
    allLLk[iter] <- oldllk


    #####################
    ## 4. Update model ##
    #####################

    # if rjmcmc
    # propose move to new model
    # accept or reject

  }
  cat("\n")
  cat("Elapsed: ", pretty_dt(difftime(Sys.time(), t0, units = "secs")), sep = "")
  cat("\n")

  return(
    c(
      inits = inits,
      priors = priors,
      allParam = allParam,
      allRates = if (allRates) allRates,
      allRateParam = if (allRateParam) allRateParam,
      allStates = allStates,
      accSwitch = if (accSwitch) allSwitch,
      accParam = accParam,
      allLen = if (allLen) allLen,
      allnLLk = allLLk,
      timing = timing
    )
  )
}

getLogPrior <- function(
    param, mu, fixPar, fixMu, priorMean, priorSD,
    rateparam, rate_priorMean, rate_priorSD, kappa, model) {
  sum(
    dnorm(
      log(param[is.na(unlist(fixPar))]),
      priorMean[1:(3 * nbStates)][is.na(unlist(fixPar))],
      priorSD[1:(3 * nbStates)][is.na(unlist(fixPar))],
      log = TRUE
    ),
    dnorm(
      log(mu[is.na(unlist(fixMu))]),
      priorMean[(3 * nbStates):(5 * nbStates)][is.na(unlist(fixMu))],
      priorSD[(3 * nbStates):(5 * nbStates)][is.na(unlist(fixMu))],
      log = TRUE
    ),
    ifelse(
      !is.na(model),
      msm::dtnorm(
        log(rateparam[1:(length(rateparam)/2)]),
        rate_priorMean,
        rate_priorSD,
        upper = log(kappa),
        log = TRUE
      ),
      0
    ),
    ifelse(
      !is.na(model),
      dunif(((length(rateparam)/2) + 1):length(rateparam), 0, 366),
      0
    )
  )
}
