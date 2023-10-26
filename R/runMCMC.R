
#' Run MCMC iterations
#'
#' @param track data.frame. Data, with columns `x`, `y`, `time`, and `ID`
#' @param nbStates integer. Number of states
#' @param nbIter integer. Number of iterations for the MCMC
#' @param inits list. Initial parameters:
#'   * `tau_pos`: vector. Initial \eqn{tau_{pos}} for each state, of length `nbStates`
#'   * `tau_vel`: vector. Initial \eqn{tau_{vel}} for each state, of length `nbStates`
#'   * `sigma`: vector. Initial \eqn{sigma} for each state, of length `nbStates` for isotropic covariance or `3 * nbStates` for anisotropic covariance
#'   * `mu`: list. Initial OUF range centre coordinate pairs `(x, y)`, with `NA` for IOU states
#'   * `Q`:
#'   * `state`: vector. Initial state sequence, length `nrow(track)`
#'   * `alpha`: vector. Length dependent on model choice
#'   * `t_alpha`: vector. Length dependent on model choice
#' @param fixed list of fixed parameters:
#'   * `tau_pos`: vector. Fixed values of \eqn{tau_{pos}} for each state with `NA` for parameters to estimate. Length `nbStates`
#'   * `tau_vel`: vector. Fixed values of \eqn{tau_{vel}} for each state with `NA` for parameters to estimate. Length `nbStates`
#'   * `sigma`: vector. Fixed values of \eqn{sigma} for each state with `NA` for parameters to estimate.
#'      Length `nbStates` for isotropic covariance, or `3 * nbStates` for anisotropic covariance
#'   * `mu`: list. Fixed OUF range centre coordinate pairs `(x, y)`, with `NA` for IOU states or pairs to estimate
#'   * `Q`: Unused
#'   * `knownStates`: vector. Known states with `NA` for those to estimate. Length `nrow(track)`
#'   * `kappa`: integer. Maximum transition rate to bound rates when they are modelled. Length 1
#' @param priors list. Parameters of prior distributions, with components:
#'   * `mean`: vector. Means for normal priors on movement parameters and Metropolis-Hastings rate parameters,
#'   of length `5 * nbStates` for movement parameters or `5 * nbStates +
#'   length(alpha)` when estimating rate parameters `alpha`, `t_alpha` by MH
#'   * `sd`: vector. Standard deviations for normal priors on movement parameters and Metropolis-Hastings rate parameters,
#'   of length `5 * nbStates` for movement parameters or `5 * nbStates +
#'   length(alpha)` when estimating rate parameters `alpha`, `t_alpha` by MH
#'   * `shape`: vector. Shapes of gamma priors for the transition rates when Gibbs sampling
#'   * `rate`: vector. Rates of gamma priors for the transition rates when Gibbs sampling
#'   * `con`: vector. Concentrations of Dirichlet priors for transition probabilities when Gibbs sampling
#' @param props list. Parameters of proposal distributions, with components:
#'   * `S`: matrix. Initial value for the lower triangular matrix of RAM algorithm, so that the covariance matrix of the proposal distribution is `SS'`.
#'   Dimensions `(5 * nbStates, 5 * nbStates)` when not modelling rate
#'   parameters (`model` is `NA`) and `(5 * nbStates + length(alpha) +
#'   length(t_alpha), 5 * nbStates + length(alpha) + length(t_alpha))`
#'   otherwise.
#'   * `updateLim`: vector. Two values: min and max length of updated state sequence
#'   * `updateProbs`: vector. Probabilities for each element of `updateLim[1]:updateLim[2]` (if `NULL`,
#'   all values are equiprobable)
#' @param tunes list. Tuning parameters, with components:
#'   * `thinStates`: integer. Thinning factor for the posterior state sequences (needed because
#'   of memory limitations)
#' @param Hmat matrix. Observation error variance (four columns, and one row for
#'   each row of data)
#' @param updateState logical. If `FALSE`, the state process is not updated (for
#'   exploratory analysis only, useful for testing single state models)
#' @param adapt integer. If `adapt` > 0, use the the Robust Adaptive Metropolis
#'   (RAM) algorithm by Vihola (2012) to update the proposal distribution for
#'   each parameter at each iteration (up to `adapt` iterations) to a target
#'   acceptance rate of 23.4%.
#' @param model experimental
#'
#' @references Michelot, T., Blackwell, P.G. (2019). State‐switching
#'   continuous‐time correlated random walks. Methods Ecol Evol, 10: 637-649.
#'   doi:10.1111/2041-210X.13154
#'
#'   Vihola, M. (2012). Robust adaptive Metropolis algorithm with coerced
#'   acceptance rate. Stat Comput, 22: 997-1008. doi:10.1007/s11222-011-9269-5
#' @examples
#' \dontrun{
#' }
#'
#' @importFrom stats dnorm runif rnorm rexp rgamma
#' @importFrom prettyunits vague_dt pretty_dt
#' @importFrom ramcmc adapt_S
#' @importFrom msm dtnorm
#' @importFrom rerddap griddap
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
  for (arg in c("tau_pos", "tau_vel")) {
    if (length(inits[[arg]]) != nbStates)
      stop("argument 'inits$", arg,"' has the wrong length, expected ", nbStates, " but got ", length(inits[[arg]]))
  }
  if (length(inits$sigma) != nbStates & length(inits$sigma) != 3 * nbStates)
    stop("argument 'inits$sigma has the wrong length, expected ", nbStates, " or ", 3 * nbStates, " but got ", length(inits$sigma))
  if (all(sapply(inits$mu, length) != rep(2, nbStates))) {
    stop("argument 'inits$mu' has the wrong length, expected ", paste(rep(2, nbStates), collapse = " "), " but got ", paste(sapply(inits$mu, length), collapse = " "))
  }
  if (length(inits$state) != nrow(track))
    stop("'inits$state' has the wrong length, expected ", nrow(track), " but got ", length(inits$state))
  if (is.null(inits$Q) & (is.null(fixed$kappa) | is.null(inits$alpha) | is.null(inits$t_alpha) | is.na(model))) {
    stop("argument 'inits$Q' is null, expected ", paste(c("fixed$kappa")[which(is.null(fixed$kappa))], c("inits$alpha", "inits$t_alpha")[which(sapply(inits[c("alpha", "t_alpha")], is.null))], c("model")[which(is.na(model))], collapse = ", "), " to be specified")
  }

  # TODO: Check Q length/dims in both init and fixed

  # Check fixed arguments and lengths
  for (arg in c("tau_pos", "tau_vel", "sigma")) {
    if (is.null(fixed[[arg]])) {
      fixed[[arg]] <- rep(NA, nbStates)
    }
    if (length(fixed[[arg]]) != length(inits[[arg]]))
      stop("argument 'fixed$", arg,"' has the wrong length, expected ", length(inits[[arg]]), " but got ", length(fixed[[arg]]))
  }
  if (is.null(fixed$mu)) {
    fixed$mu <- rep(list(c(NA,NA)), nbStates)
  } else {
    if (all(sapply(fixed$mu, length) != rep(2, nbStates))) {
      stop("argument 'fixed$mu' has the wrong length, expected ", paste(rep(2, nbStates), collapse = " "), " but got ", paste(sapply(fixed$mu, length), collapse = " "))
    }
  }
  if (is.null(fixed$knownStates) | is.na(fixed$knownStates)) {
    fixed$knownStates <- rep(NA, nrow(track))
  } else {
    if (length(fixed$knownStates) != nrow(track)) {
      stop("argument 'fixed$knownStates' has the wrong length, expected ", nrow(track), " but got ", length(fixed$knownStates))
    }
  }
  if (length(fixed$kappa) != 1)
    stop("argument 'fixed$kappa' has the wrong length, expected ", 1, " but got ", length(fixed$kappa))

  nbParam <- 4 + length(inits$sigma) / nbStates

  # Check priors arguments and lengths
  # TODO smarter calculation when args are fixed
  if (is.null(priors$mean) | is.null(priors$sd))
    stop("argument 'priors' missing: ", paste(c("mean", "sd")[which(sapply(priors[c("mean", "sd")], is.null))], collapse = ", "))
  for (arg in c("mean", "sd")) {
    # We dont check for t_alpha here because it has a uniform prior specified elsewhere
    if (is.na(model) & length(priors[[arg]]) != nbParam * nbStates)
      stop("argument 'priors$", arg,"' has the wrong length, expected ", nbParam * nbStates, " but got ", length(priors[[arg]]))
    if (!is.na(model) & length(priors[[arg]]) != nbParam * nbStates + length(inits$alpha))
      stop("argument 'priors$", arg,"' has the wrong length, expected ", nbParam * nbStates + length(inits$alpha), " but got ", length(priors[[arg]]))
  }
  if (is.na(model) & (is.null(priors$shape) | is.null(priors$rate) | is.null(priors$con)))
    stop("argument 'model' is NA, but argument 'priors' missing: ", paste(c("shape", "rate", "con")[which(sapply(priors[c("shape", "rate", "con")], is.null))], collapse = ", "))
  if (is.na(model)) {
    for (arg in c("shape", "rate", "con")) {
      if (length(priors[[arg]]) != nbStates)
        stop("argument 'priors$", arg,"' has the wrong length, expected ", nbStates, " but got ", length(priors[[arg]]))
    }
  }

  # Check props arguments and lengths
  if (is.null(props$S) | (updateState & is.null(props$updateLim)))
    stop("argument 'props' missing: ", paste(c("S")[which(is.null(props$S))], c("S")[which(updateState & is.null(props$S))], collapse = ", "))
  if (is.na(model) & all(dim(props$S) != c(nbParam * nbStates, nbParam * nbStates)))
    stop("argument 'props$S' has the wrong dimensions, expected ", paste(c(nbParam * nbStates, nbParam * nbStates), collapse = ", "), " but got ", paste(dim(props$S), collapse = ", "))
  if (!is.na(model) & all(dim(props$S) != c(nbParam * nbStates + length(inits$alpha) + length(inits$t_alpha), nbParam * nbStates + length(inits$alpha) + length(inits$t_alpha))))
    stop("argument 'props$S' has the wrong dimensions, expected ", paste(c(nbParam * nbStates + length(inits$alpha) + length(inits$t_alpha), nbParam * nbStates + length(inits$alpha) + length(inits$t_alpha)), collpase = ", "), " but got ", paste(dim(props$S), collapse = ", "))
  if (updateState & "list" %in% class(props$updateLim) & all(sapply(props$updateLim, length) != rep(2, length(unique(track$ID)))))
    stop("argument 'props$updateLim' has the wrong length, expected ", paste(rep(2, length(unique(track$ID))), collapse = " "), " but got ", paste(sapply(props$updateLim, length), collapse = " "))
  if (updateState & !"list" %in% class(props$updateLim) & length(props$updateLim) != 2)
    stop("argument 'props$updateLim' has the wrong length, expected ", 2, " but got ", length(props$updateLim))
  if (updateState & !is.null(props$updateProbs) & class(props$updateLim) != class(props$updateProbs))
    stop("argument 'props$updateProbs' provided but of wrong type, expected ", class(props$updateLim), " but got ", class(props$updateProbs))
  if (updateState & !is.null(props$updateProbs) & "list" %in% class(props$updateProbs) & all(sapply(props$updateLim, length) != sapply(props$updateProbs, length)))
    stop("argument 'props$updateProbs' has the wrong length, expected ", paste(sapply(props$updateLim, length), collapse = " "), " but got ", paste(sapply(props$updateProbs, length), collapse = " "))
  if (updateState & !is.null(props$updateProbs) & !"list" %in% class(props$updateProbs) & length(props$updateProbs) != 2)
    stop("argument 'props$updateProbs' has the wrong length, expected ", 2, " but got ", length(props$updateProbs))

  # TODO: Use missing(model) instead?
  if (!updateState & (!is.null(inits$Q) | !is.null(fixed$Q)))
    warning("argument 'updateState' is FALSE, ignoring 'inits$Q', 'fixed$Q', 'priors$shape', 'priors$rate', 'priors$con'")
  if (!is.na(model) & (!is.null(inits$Q)))
     warning("argument 'model' is not NA, ignoring 'inits$Q'")
  if (!is.na(model) & (!is.null(priors$shape) | !is.null(priors$rate) | !is.null(priors$con)))
    warning("argument 'model' is not NA, ignoring ", paste(c("priors$shape", "priors$rate", "priors$con")[which(sapply(priors[c("shape", "rate", "con")], is.null))], collapse = ", "))
  if (is.na(model) & all(dim(props$S) > c(nbParam * nbStates, nbParam * nbStates)))
    warning("argument 'model' is NA, ignoring extra dimensions of 'props$S")

  # TODO check any NA in Hmat
  # TODO abstract out model (accept function as argument for model and priors)

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
  priorMean <- priors$mean[1:(nbParam * nbStates)]
  priorSD <- priors$sd[1:(nbParam * nbStates)]
  priorShape <- priors$shape
  priorRate <- priors$rate
  priorCon <- priors$con

  # check if rate matrix provided
  # TODO just check model here?
  if (!is.null(inits$Q) & length(inits$Q) == length(unique(track$ID)) & "list" %in% class(inits$Q)) {
    Q <- inits$Q
    names(Q) <- unique(track$ID)
    kappa <- fixed$kappa
    ratePriorMean <- NULL
    ratePriorSD <- NULL
    rateparam <- NULL
  } else { # if no rate matrix provided, rate params must be
    Q <- NULL
    kappa <- fixed$kappa
    rateS <- props$S[(nbParam * nbStates + 1):nrow(props$S), (nbParam * nbStates + 1):ncol(props$S)]
    ratePriorMean <- priors$mean[(nbParam * nbStates + 1):length(priors$mean)]
    ratePriorSD <- priors$sd[(nbParam * nbStates + 1):length(priors$sd)]
    rateparam <- c(inits$alpha, inits$t_alpha)
  }

  # unpack proposal parameters
  S <- props$S[1:(nbParam * nbStates), 1:(nbParam * nbStates)]
  # TODO: CLEAN THIS UP?
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
      rep(1, length(updateLim[[i]][1]:updateLim[[i]][2])) / length(updateLim[[i]][1]:updateLim[[i]][2])
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

  names(updateLim) <- ids
  if (exists("updateProbs")) names(updateProbs) <- ids

  for (id in ids) {
    nbObs <- nrow(track[which(track$ID == id), ])
    obs[[id]] <- cbind(
      "x" = track[which(track$ID == id), "x"],
      "y" = track[which(track$ID == id), "y"],
      "time" = track[which(track$ID == id), "time"],
      "ID" = as.numeric(track[which(track$ID == id), "ID"]),
      "state" = state[which(track$ID == id)],
      "group" = ifelse("group" %in% colnames(track), track[which(track$ID == id), "group"], 1)
    )
    known[[id]] <- knownStates[which(track$ID == id)]
    #colnames(obs[[id]]) <- c("x", "y", "time", "ID",  "state")
    indSwitch <- which(obs[[id]][-1, "state"] != obs[[id]][-nbObs, "state"]) + 1
    switch[[id]] <- cbind(
      "time" = obs[[id]][indSwitch, "time"] - 0.001,
      "state" = rle(obs[[id]][ , "state"])$values[-1])
    #colnames(switch[[id]]) <- c("time", "state")

    if (!all(is.na(switch[[id]]))) {
      data.list[[id]] <- rbind(
        obs[[id]][ , c( "x", "y", "time", "ID",  "state", "group")],
        cbind(
          "x" = NA,
          "y" = NA,
          "time" = switch[[id]][ , "time"],
          "ID" = id,
          "state" = switch[[id]][ , "state"],
          "group" = ifelse("group" %in% colnames(track), track[which(track$ID == id), "group"], 1)
        )
      )
      data.list[[id]] <- data.list[[id]][order(data.list[[id]][ , "time"]), ]
    } else {
      data.list[[id]] <- obs[[id]][ , c("x", "y", "time", "ID",  "state", "group")]
    }
  }

  # flatten data
  data.mat <- do.call("rbind", data.list)
  data.mat <- data.mat[ , c("x", "y", "time", "ID", "state", "group")]

  # initialise Hmat (rows of 0s for transitions)
  HmatAll <- matrix(0, nrow(data.mat), 4)
  HmatAll[which(!is.na(data.mat[ , "x"])), ] <- Hmat

  # initial likelihood
  oldllk <- kalman_rcpp(data = data.mat, nbStates = nbStates, param = param, fixmu = mu, Hmat = HmatAll)$llk

  # initial log-prior
  oldlogprior <- getLogPrior(
    param, mu, fixPar, fixMu, priorMean, priorSD,
    rateparam, ratePriorMean, ratePriorSD, kappa, model
  )

  ###############################
  ## Loop over MCMC iterations ##
  ###############################

  allParam <- matrix(NA, nrow = nbIter, ncol = nbParam * nbStates)
  colnames(allParam) <- c(
    paste("tau_pos[", 1:nbStates, "]", sep = ""),
    paste("tau_vel[", 1:nbStates, "]", sep = ""),
    switch((nbParam == 7) + 1, paste("sigma[", 1:nbStates, "]", sep = ""), paste(c("sigma_x[", "sigma_y[", "sigma_xy["), rep(1:nbStates, each = 3), rep("]", nbStates * 3), sep = "")),
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

    } else {
      allRateParam <- matrix(NA, nrow = nbIter, ncol = length(rateparam))
    }
  }


  timing <- matrix(NA, nrow = nbIter / thinStates, ncol = 2)

  t0 <- Sys.time()
  for (iter in 1:nbIter) {
    if (iter < 100 | iter %% 100 == 0) {
      cat(
        "\33[2K\rIteration ", iter, "/", nbIter, "... ",
        vague_dt(difftime(Sys.time(), t0, units = "secs") / iter * (nbIter - iter), "short"),
        " remaining (est)",
        ifelse(
          exists("accSwitch"),
          paste0(" -- accSwitch = ", round(sum(accSwitch) / iter * 100), "%"),
          ""
        ),
        " -- accParam = ", round(sum(accParam) / iter * 100), "%", sep = ""
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
        newRateParams <- proposeParams(rateparam, NA, rateS)
      } else {
        newRateParams <- NULL
      }

      logHR <- -Inf

      try({ # catch autoreject signal
        for (id in ids) {
          #TODO: FOR GHOST TOWN, WE NEED ALL DATA, NOT JUST OBS
          upState <- updateState(
            obs = obs[[id]],
            nbStates = nbStates,
            knownStates = known[[id]],
            switch = switch[[id]],
            updateLim = updateLim[[id]],
            param = param,
            mu = unlist(mu),
            Hmat = HmatAll[which(data.mat[,"ID"] == id),],
            updateProbs = updateProbs[[id]],
            Q = switch(is.null(Q) + 1, Q[[id]], NULL), #https://www.r-bloggers.com/2017/02/use-switch-instead-of-ifelse-to-return-a-null/
            rateparam = newRateParams[[2]],
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
        newllk <- kalman_rcpp(data = newData.mat, nbStates = nbStates, param = param, fixmu = unlist(mu), Hmat = newHmatAll)$llk
        newlogprior <- getLogPrior(
          param, mu, fixPar, fixMu, priorMean, priorSD,
          newRateParams[[2]], ratePriorMean, ratePriorSD, kappa, model
        )
        logHR <- newllk + newlogprior - oldllk - oldlogprior

        if (log(runif(1)) < logHR) {
          # Accept new state sequence
          accSwitch[iter] <- 1
          switch <- newSwitch
          data.list <- newData.list
          obs <- lapply(data.list, function(data) { data[!is.na(data[ , "x"]), ] })
          oldllk <- newllk
          rateparam <- newRateParams[[2]]
          oldlogprior <- newlogprior
          HmatAll <- newHmatAll
        }

        #if (adapt & iter >= 1000 & iter <= adapt) {
        if (!is.na(model) & adapt & iter > 1 & iter <= adapt) {
          rateS <- adapt_S(rateS, newRateParams[[1]], min(1, exp(logHR)), iter)
        }
      #}, silent = TRUE)
      })

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
    }

    ###################################
    ## 2. Update movement parameters ##
    ###################################

    pass = F
    while (!pass) {  # ensure tau_p >= tau_v

      # On working scale [-Inf,Inf]
      newParams <- proposeParams(param, fixPar, S[1:length(param), 1:length(param)])
      newMu <- proposeParams(mu, fixMu, S[(length(param) + 1):nrow(S), (length(param) + 1):ncol(S)])
      # NB we could bound mu to -180,180 -90,90 with a different dist in proposeParams() but would need projected bounds

      # hack to ensure tau_pos >= tau_vel (does actually limit models we can test)
      # is there potential to get stuck here?
      if (all(newParams[[2]][1:nbStates] >= newParams[[2]][(nbStates + 1):(2 * nbStates)])) {
        pass <- T
      }
    }

    # Calculate acceptance ratio
    data.mat <- do.call("rbind", data.list)
    data.mat <- data.mat[ , c("x", "y", "time", "ID", "state")]
    newlogprior <- getLogPrior(
      newParams[[2]], newMu[[2]], fixPar, fixMu, priorMean, priorSD,
      rateparam, ratePriorMean, ratePriorSD, kappa, model
    )
    kalman <- kalman_rcpp(data = data.mat, nbStates = nbStates, param = newParams[[2]], fixmu = newMu[[2]], Hmat = HmatAll)
    newllk <- kalman$llk
    #mu <- as.vector(t(kalman$mu))
    logHR <- newllk + newlogprior - oldllk - oldlogprior

    if (log(runif(1)) < logHR) {
      # Accept new parameter values
      accParam[iter] <- 1
      param <- newParams[[2]]
      mu <- newMu[[2]]
      oldllk <- newllk
      oldlogprior <- newlogprior
    }

    #if (adapt & iter >= 1000 & iter <= adapt) {
    if (adapt & iter <= adapt) {
      #S[is.na(unlist(fixPar)), is.na(unlist(fixPar))] <- adapt_S(S[is.na(unlist(fixPar)), is.na(unlist(fixPar))], param_u[is.na(unlist(fixPar))], min(1, exp(logHR)), iter)
      # calculate S by state instead
      for (i in 1:nbStates) {
        #TODO this will change with sigma = 3
        paramindex <- seq(i, length(param), by = nbStates)
        S[paramindex, paramindex][is.na(unlist(fixPar)[paramindex]), is.na(unlist(fixPar)[paramindex])] <- adapt_S(S[paramindex,paramindex][is.na(unlist(fixPar)[paramindex]), is.na(unlist(fixPar)[paramindex])], newParams[[1]][paramindex][is.na(unlist(fixPar)[paramindex])], min(1, exp(logHR)), iter)
        muindex <- c(i * 2 - 1, i * 2)
        if (any(is.na(unlist(fixMu)[muindex])))
          S[length(param) + muindex, length(param) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])] <- adapt_S(S[length(param) + muindex, length(param) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])], newMu[[1]][muindex][is.na(unlist(fixMu)[muindex])], min(1, exp(logHR)), iter)
      }
    }

    #########################
    ## Save posterior draw ##
    #########################
    allParam[iter,] <- cbind(matrix(param, ncol = length(param)), matrix(mu, ncol =  2 * nbStates))
    if (updateState & !is.null(Q)) {
      allRates[iter, , ] <- matrix(unlist(lapply(Q, function(q){ q[!diag(nbStates)]})), ncol = length(ids), nrow = nbStates * (nbStates - 1))
    } else if (updateState) {
      allRateParam[iter, ] <- rateparam
    }

    if (iter %% thinStates == 0) {
      if (updateState)
        allStates[iter / thinStates, ] <- unlist(lapply(obs, function(ob) { ob[ , "state"] }))
      timing[iter / thinStates, ] <- c(iter, Sys.time())
    }
    allLLk[iter] <- oldllk


    #####################
    ## 4. Update model ##
    #####################

    # if rjmcmc
    # propose move to new model
      #if nMu = 1
        # split
      # else
        # split or merge (random bernoulli)
        # if split
          # where to split?
        # if merge
          # where to merge?
    # accept or reject
    # NB split or merge update mu, Q, params, posssibly rateMean, rateSD, S...
    # https://www.gen.dev/tutorials/rj/tutorial#split-merge

  }
  cat("\n")
  cat("Elapsed: ", pretty_dt(difftime(Sys.time(), t0, units = "secs")), sep = "")
  cat("\n")

  return(
    c(
      list(inits = inits),
      list(priors = priors),
      list(allParam = allParam),
      if (exists("allRates")) list(allRates = allRates),
      if (exists("allRateParam")) list(allRateParam = allRateParam),
      if (exists("allStates")) list(allStates = allStates),
      if (exists("accSwitch")) list(accSwitch = accSwitch),
      list(accParam = accParam),
      if (exists("allLen")) list(allLen = allLen),
      list(allnLLk = allLLk),
      list(timing = timing)
    )
  )
}

proposeParams <- function(param, fixedParams, S) {
  sign <- ifelse(param < 0, -1, 1)

  # On working scale [-Inf, Inf]
  u <- rnorm(length(param))
  thetas <- suppressWarnings(c(log(abs(param))) + as.vector(S %*% u)) #In log(param) : NaNs produced

  if (all(is.na(fixedParams))) {
    fixedParams <- sapply(param, function(x) {rep(NA, length(x))})
  }

  # On natural scale [0, Inf]
  thetasprime <- unlist(fixedParams)
  thetasprime[is.na(unlist(fixedParams))] <- exp(thetas[is.na(unlist(fixedParams))]) * sign[is.na(unlist(fixedParams))]
  return(list(u, thetasprime))
}

getLogPrior <- function(
    param, mu, fixPar, fixMu, priorMean, priorSD,
    rateparam, ratePriorMean, ratePriorSD, kappa, model) {
  return(
    sum(
      dnorm(
        log(param[is.na(unlist(fixPar))]),
        priorMean[1:length(unlist(fixPar))][is.na(unlist(fixPar))],
        priorSD[1:length(unlist(fixPar))][is.na(unlist(fixPar))],
        log = TRUE
      ),
      dnorm( # mu are not on log scale (can be negative)
        mu[is.na(unlist(fixMu)) & !is.na(mu)],
        priorMean[(length(unlist(fixPar)) + 1):length(priorMean)][is.na(unlist(fixMu)) & !is.na(mu)],
        priorSD[(length(unlist(fixPar)) + 1):length(priorSD)][is.na(unlist(fixMu)) & !is.na(mu)],
        log = TRUE
      ),
      ifelse(
        !is.na(model),
        msm::dtnorm(
          log(rateparam[1:(length(rateparam)/2)]),
          ratePriorMean,
          ratePriorSD,
          upper = log(kappa),
          log = TRUE
        ),
        0
      ),
      ifelse(
        !is.na(model),
        dunif(
          log(rateparam[((length(rateparam)/2) + 1):length(rateparam)]),
          0,
          40, # NB: CHANGE THIS!
          log = TRUE
        ),
        0
      )
    )
  )
}
