#' Run MCMC iterations
#'
#' @param track data.frame. Data, with columns `x`, `y`, `time`, and `ID`
#' @param nbStates integer. Number of states
#' @param nbIter integer. Number of iterations for the MCMC
#' @param inits list. Initial parameters:
#'   * `tau_pos`: vector. Initial \eqn{tau_{pos}} for each state,
#'      of length `nbStates`
#'   * `tau_pos_sd`: vector (optional). If provided, initial \eqn{sigma}
#'      hyperparameters for individual random effects on tau_pos for each state,
#'      tau_pos argument will be used as \eqn{mu} hyperparameters,
#'      of length `nbStates`
#'   * `tau_vel`: vector. Initial \eqn{tau_{vel}} for each state,
#'      of length `nbStates`
#'   * `tau_vel_sd`: vector (optional). If provided, initial \eqn{sigma}
#'      hyperparameters for individual random effects on tau_vel for each state,
#'      tau_vel argument will be used as \eqn{mu} hyperparameters,
#'      of length `nbStates`
#'   * `sigma`: vector. Initial \eqn{sigma} for each state,
#'      of length `nbStates` for isotropic covariance
#'      or `3 * nbStates` for anisotropic covariance
#'   * `sigma_sd`: vector (optional). If provided, initial \eqn{sigma}
#'      hyperparameters for individual random effects on sigma for each state,
#'      sigma argument will be used as \eqn{mu} hyperparameters,
#'      of length `nbStates` or `3 * nbStates`
#'   * `mu`: list. Initial range centre coordinate pairs `(x, y)`,
#'      with `NA` for OU states
#'   * `Q`:
#'   * `state`: vector. Initial state sequence, length `nrow(track)`
#'   * `rateparam`: vector. Length dependent on model choice
#' @param fixed list of fixed parameters:
#'   * `tau_pos`: vector. Fixed values of \eqn{tau_{pos}} for each state
#'      with `NA` for parameters to estimate. Length `nbStates`
#'   * `tau_vel`: vector. Fixed values of \eqn{tau_{vel}} for each state
#'      with `NA` for parameters to estimate. Length `nbStates`
#'   * `sigma`: vector. Fixed values of \eqn{sigma} for each state
#'      with `NA` for parameters to estimate.
#'      Length `nbStates` for isotropic covariance,
#'      or `3 * nbStates` for anisotropic covariance
#'   * `mu`: list. Fixed range centre coordinate pairs `(x, y)`,
#'      with `NA` for IOU states or pairs to estimate
#'   * `Q`: Unused
#'   * `knownStates`: vector. Known states with `NA` for those to estimate.
#'      Length `nrow(track)`
#'   * `kappa`: integer. Maximum transition rate to bound rates
#'      when they are modelled. Length 1
#' @param priors list. Parameters of prior distributions, with components:
#'   * `func`: list of distribution function names,
#'      of length `5 * nbStates`
#'      or `7 * nbStates` for movement parameters + `length(rateparam)`
#'      when estimating rate parameters `rateparam` by MH
#'   * `args`: list of lists containing function args for distribution
#'      functions, of length `5 * nbStates`
#'      or `7 * nbStates` for movement parameters + `length(rateparam)`
#'      when estimating rate parameters `rateparam` by MH
#'   * `shape`: vector. Shapes of gamma priors for the transition rates
#'      when Gibbs sampling
#'   * `rate`: vector. Rates of gamma priors for the transition rates
#'      when Gibbs sampling
#'   * `con`: vector. Concentrations of Dirichlet priors for
#'      transition probabilities when Gibbs sampling
#' @param props list. Parameters of proposal distributions, with components:
#'   * `S`: matrix. Initial value for the lower triangular matrix of
#'      RAM algorithm, so that the covariance matrix of the
#'      proposal distribution is `SS'`.
#'      Dimensions `(5 * nbStates, 5 * nbStates)` when not modelling rate
#'      parameters (`model` is `NA`) and
#'      `(5 * nbStates + length(rateparam), 5 * nbStates + length(rateparam))`
#'      otherwise.
#'   * `updateLim`: vector. Two values: min and max length of
#'      updated state sequence
#'   * `updateProbs`: vector. Probabilities for each element
#'      of `updateLim[1]:updateLim[2]` (if `NULL`, all values are equiprobable)
#' @param tunes list. Tuning parameters, with components:
#'   * `thinStates`: integer. Thinning factor for the posterior state sequences
#'     (useful for reducing memory)
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
#'
#' }
#'
#' @importFrom stats dnorm runif rnorm rexp rgamma
#' @importFrom msm dtnorm
#' @importFrom prettyunits vague_dt pretty_dt
#' @importFrom ramcmc adapt_S
#' @importFrom rerddap griddap
#' @export
#'
#' @useDynLib MSctmm
runMCMC <- function(
    track, nbStates, nbIter, burnin, inits, fixed, priors,
    props, tunes, Hmat, updateState = TRUE, adapt = FALSE, model = NA,
    debug = FALSE) {
  options(warn = 1) # show warnings as they occur

  # Check track df
  track |> stopIfNotDataFrame()
  track |> stopIfMissingColumns(c("x", "y", "time"))

  if (isMissing(track, "ID")) {
    warning("argument 'track' should have column 'ID', assuming data is one individual")
    track$ID <- 1
  } else {
    track$ID <- as.numeric(as.factor(track$ID))
  }

  # Check inits arguments and lengths
  inits |> stopIfMissingObjects(c("tau_pos", "tau_vel", "sigma", "mu", "state"))
  inits |> stopIfLengthNotEquals(c("tau_pos", "tau_vel"), nbStates)
  inits |> stopIfLengthNotEquals("sigma", c(nbStates, 3 * nbStates))
  inits |> stopIfNestedLengthNotEquals("mu", rep(2, nbStates))
  inits |> stopIfLengthNotEquals("state", nrow(track))

  if (all(updateState, isMissing(inits, "Q"), any(
    isMissing(fixed, "kappa"),
    isMissing(inits, "rateparam"), is.na(model)
  ))) {
    stop(
      "argument 'inits$Q' is null, expected ",
      paste(
        c("fixed$kappa")[which(is.null(fixed$kappa))],
        c("inits$rateparam")[which(sapply(inits[c("rateparam")], is.null))],
        c("model")[which(is.na(model))],
        collapse = ", "
      ),
      " to be specified"
    )
  }

  # Check for random effects hyperparameters
  if (all(!isMissing(inits, c("tau_pos_sd", "tau_vel_sd", "sigma_sd")))) {
    # TODO: check priors
    randomEffects <- TRUE
  } else {
    randomEffects <- FALSE
  }

  # TODO: Check Q length/dims in both init and fixed

  # Check fixed arguments and lengths
  for (arg in c("tau_pos", "tau_vel", "sigma")) {
    if (isMissing(fixed, arg)) {
      fixed[[arg]] <- rep(NA, nbStates)
    }
    fixed |> stopIfLengthNotEquals(arg, length(inits[[arg]]))
  }
  if (isMissing(fixed, "mu")) {
    fixed$mu <- rep(list(c(NA, NA)), nbStates)
  } else {
    fixed |> stopIfNestedLengthNotEquals("mu", rep(2, nbStates))
  }
  if (isMissing(fixed, "knownStates") || is.na(fixed$knownStates)) {
    fixed$knownStates <- rep(NA, nrow(track))
  } else {
    fixed |> stopIfLengthNotEquals("knownStates", nrow(track))
  }
  fixed |> stopIfLengthNotEquals("kappa", 1)

  # number movement parameters (including mu) per state
  nbParam <- 4 + length(inits$sigma) / nbStates
  # number of hyperparams (should be 2x movement params, excluding mu) if estimating random effects
  nbHyperParam <- ifelse(randomEffects, (2 + length(inits$sigma) / nbStates), 0)

  # Check priors arguments and lengths
  # TODO smarter calculation when args are fixed
  priors |> stopIfMissingObjects(c("func", "args"))
  if (is.na(model)) {
    if (!randomEffects) {
      priors |> stopIfLengthNotEquals(c("func", "args"), nbParam * nbStates)
    } else {
      priors |> stopIfLengthNotEquals(c("func", "args"), nbParam * nbStates + nbHyperParam * nbStates)
    }
  } else {
    if (!randomEffects) {
      priors |>
        stopIfLengthNotEquals(
          c("func", "args"),
          nbParam * nbStates + length(inits$rateparam)
        )
    } else {
      priors |>
        stopIfLengthNotEquals(
          c("func", "args"),
          nbParam * nbStates + length(inits$rateparam) + nbHyperParam * nbStates
        )
    }
  }

  if (is.na(model) &&
    any(
      isMissing(priors, "shape"),
      isMissing(priors, "rate"),
      isMissing(priors, "con")
    )) {
    stop(
      "argument 'model' is NA, but argument 'priors' missing: ",
      paste(c("shape", "rate", "con")[which(sapply(priors[c("shape", "rate", "con")], is.null))], collapse = ", ")
    )
  }
  if (is.na(model)) {
    priors |> stopIfLengthNotEquals(c("shape", "rate", "con"), nbStates)
  }

  # TODO: check that all inits names exist in priors$funcs and priors$args

  # Check props arguments and lengths
  if (updateState) props |> stopIfMissingObjects("updateLim")

  props |> stopIfMissingObjects("S")
  # TODO: handle dimension when randomEffects
  # should have an S triangle for hyperparams for each state

  if (is.na(model)) {
    props |> stopIfDimNotEquals("S", rep(nbParam * nbStates + nbHyperParam * nbStates, 2))
  } else {
    props |>
      stopIfDimNotEquals(
        "S",
        rep(nbParam * nbStates + nbHyperParam * nbStates + length(inits$rateparam), 2)
      )
  }
  if (all(updateState, "list" %in% class(props$updateLim))) {
    props |>
      stopIfNestedLengthNotEquals("updateLim", rep(2, length(unique(track$ID))))
  }
  if (all(updateState, !"list" %in% class(props$updateLim))) {
    props |> stopIfLengthNotEquals("updateLim", 2)
  }
  if (all(
    updateState, !is.null(props$updateProbs),
    class(props$updateLim) != class(props$updateProbs)
  )) {
    stop(
      "argument 'props$updateProbs' provided but of wrong type, expected ",
      class(props$updateLim),
      " but got ",
      class(props$updateProbs)
    )
  }
  if (all(
    updateState, !isMissing(props, "updateProbs"),
    "list" %in% class(props$updateProbs)
  )) {
    props |>
      stopIfNestedLengthNotEquals(
        "updateLim",
        sapply(props$updateProbs, length)
      )
  }
  if (all(
    updateState, !isMissing(props, "updateProbs"),
    !"list" %in% class(props$updateProbs)
  )) {
    props |> stopIfLengthNotEquals("updateProbs", 2)
  }

  # TODO: Use missing(model) instead?
  if (!updateState && any(!isMissing(inits, "Q"), !isMissing(fixed, "Q"))) {
    warning(
      "argument 'updateState' is FALSE, ignoring 'inits$Q', 'fixed$Q', 'priors$shape', 'priors$rate', 'priors$con'"
    )
  }
  if (!is.na(model) && !isMissing(inits, "Q")) {
    warning("argument 'model' is not NA, ignoring 'inits$Q'")
  }
  if (!is.na(model) &&
    any(!isMissing(priors, c("shape", "rate", "con")))) {
    warning(
      "argument 'model' is not NA, ignoring ",
      paste(c(
        "priors$shape", "priors$rate", "priors$con"
      )[which(!sapply(priors[c("shape", "rate", "con")], is.null))], collapse = ", ")
    )
  }
  if (is.na(model) &&
    all(dim(props$S) > c(nbParam * nbStates, nbParam * nbStates))) { # TODO: update this to reflect hyperparams if random effects
    warning("argument 'model' is NA, ignoring extra dimensions of 'props$S")
  }

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
  if (randomEffects) {
    # Need full set of params for each individual
    hyperparam_mu <- unlist(inits[c("tau_pos", "tau_vel", "sigma")])
    hyperparam_sd <- unlist(inits[c("tau_pos_sd", "tau_vel_sd", "sigma_sd")])
    param <- lapply(unique(track$ID), function(id) {
      # TODO: hide warning here for Inf?
      p <- rtnorm(length(hyperparam_mu), mean = hyperparam_mu, sd = hyperparam_sd, lower = sapply(names(hyperparam_mu), function(par) {
        priors$args[[par]]$min
      }), upper = sapply(names(hyperparam_mu), function(par) {
        priors$args[[par]]$max
      }))
      names(p) <- names(hyperparam_mu)
      # overwrite fixed params
      p[names(which(!is.na(unlist(fixPar))))] <- unlist(fixPar)[names(which(!is.na(unlist(fixPar))))]
      return(p)
    })
  } else {
    param <- unlist(inits[c("tau_pos", "tau_vel", "sigma")])
    # overwrite fixed params
    param[names(which(!is.na(unlist(fixPar))))] <- unlist(fixPar)[names(which(!is.na(unlist(fixPar))))]
  }

  mu <- unlist(inits$mu)
  state <- inits$state

  # overwrite other fixed params
  mu[which(!is.na(unlist(fixMu)))] <-
    unlist(fixMu)[which(!is.na(unlist(fixMu)))]
  state[which(!is.na(knownStates))] <-
    knownStates[which(!is.na(knownStates))]

  # ensure tau_p > tau_v
  # TODO subset by name not index
  if (randomEffects) {
    for (id in unique(track$ID)) {
      param[[id]][1:(2 * nbStates)] <- c(
        pmax(param[[id]][1:nbStates], param[[id]][(nbStates + 1):(2 * nbStates)]),
        pmin(param[[id]][1:nbStates], param[[id]][(nbStates + 1):(2 * nbStates)])
      )
    }
  } else {
    param[1:(2 * nbStates)] <- c(
      pmax(param[1:nbStates], param[(nbStates + 1):(2 * nbStates)]),
      pmin(param[1:nbStates], param[(nbStates + 1):(2 * nbStates)])
    )
  }

  # unpack prior parameters
  # priorFunc <- priors$func[1:(nbParam * nbStates)]
  # priorArgs <- priors$args[1:(nbParam * nbStates)]
  priorFunc <- priors$func[!str_starts(names(priors$func), "rateparam")]
  priorArgs <- priors$args[!str_starts(names(priors$args), "rateparam")]
  priorShape <- priors$shape
  priorRate <- priors$rate
  priorCon <- priors$con

  # check if rate matrix provided
  # TODO just check model here?
  if (!updateState) {
    ratePriorFunc <- NULL
    ratePriorArgs <- NULL
    rateparam <- NULL
  } else if (all(
    is.na(model), !isMissing(inits, "Q"),
    length(inits$Q) == length(unique(track$ID)),
    "list" %in% class(inits$Q)
  )) {
    Q <- inits$Q
    names(Q) <- unique(track$ID)
    kappa <- fixed$kappa
    ratePriorFunc <- NULL
    ratePriorArgs <- NULL
    rateparam <- NULL
  } else {
    # if no rate matrix provided, rate params must be
    Q <- NULL
    kappa <- fixed$kappa
    rateS <-
      props$S[(nbParam * nbStates + nbHyperParam * nbStates + 1):nrow(props$S), (nbParam * nbStates + nbHyperParam * nbStates + 1):ncol(props$S)]
    # ratePriorFunc <-
    #  priors$func[(nbParam * nbStates + 1):length(priors$func)]
    # ratePriorArgs <-
    #  priors$args[(nbParam * nbStates + 1):length(priors$args)]
    ratePriorFunc <- priors$func[str_starts(names(priors$func), "rateparam")]
    ratePriorArgs <- priors$args[str_starts(names(priors$args), "rateparam")]
    rateparam <- inits$rateparam
  }

  # unpack proposal parameters
  S <- props$S[1:(nbParam * nbStates + nbHyperParam * nbStates), 1:(nbParam * nbStates + nbHyperParam * nbStates)]

  # initialize random effect scale matrix from props
  if (randomEffects) {
    randomS <- lapply(unique(track$ID), function(id) {
      S[1:((nbParam - 2) * nbStates), 1:((nbParam - 2) * nbStates)]
    })
  }

  # TODO: CLEAN THIS UP?
  if ("list" %in% class(props$updateLim) &&
    length(props$updateLim) == length(unique(track$ID))) {
    updateLim <- lapply(seq_along(unique(track$ID)), function(i) {
      lim <-
        ceiling(props$updateLim[[i]] * nrow(track[which(track$ID == unique(track$ID)[i]), ]))
      if (lim[1] <= 3) {
        lim <- lim + 1
      }
      if (lim[1] == lim[2]) {
        lim[2] <- lim[2] + 1
      }
      lim
    })
  } else {
    updateLim <- lapply(unique(track$ID), function(id) {
      lim <-
        ceiling(props$updateLim * nrow(track[which(track$ID == id), ]))
      if (lim[1] <= 3) {
        lim <- lim + 1
      }
      if (lim[1] == lim[2]) {
        lim[2] <- lim[2] + 1
      }
      lim
    })
  }
  if (!is.null(props$updateProbs) &&
    "list" %in% class(props$updateProbs) &&
    length(props$updateProbs) == length(unique(track$ID))) {
    updateProbs <- props$updateProbs
    # } else {
    #  updateProbs <- rep(list(props$updateProbs), length(unique(track$ID)))
  } else if (is.null(props$updateProbs)) {
    updateProbs <- lapply(seq_along(unique(track$ID)), function(i) {
      rep(1, length(updateLim[[i]][1]:updateLim[[i]][2])) / length(updateLim[[i]][1]:updateLim[[i]][2])
    })
  }

  # unpack tuning parameters
  thinStates <-
    ifelse(is.null(tunes$thinStates), 1, tunes$thinStates)

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
  if (exists("updateProbs")) {
    names(updateProbs) <- ids
  }

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
    # colnames(obs[[id]]) <- c("x", "y", "time", "ID",  "state")
    indSwitch <-
      which(obs[[id]][-1, "state"] != obs[[id]][-nbObs, "state"]) + 1
    switch[[id]
    ] <- cbind(
      "time" = obs[[id]][indSwitch, "time"] - 0.001,
      "state" = rle(obs[[id]][, "state"])$values[-1]
    )
    # colnames(switch[[id]]) <- c("time", "state")

    if (!all(is.na(switch[[id]
    ]))) {
      data.list[[id]] <- rbind(
        obs[[id]][, c("x", "y", "time", "ID", "state", "group")],
        cbind(
          "x" = NA,
          "y" = NA,
          "time" = switch[[id]
          ][, "time"],
          "ID" = id,
          "state" = switch[[id]
          ][, "state"],
          "group" = ifelse("group" %in% colnames(track), track[which(track$ID == id), "group"], 1)
        )
      )
      data.list[[id]] <-
        data.list[[id]][order(data.list[[id]][, "time"]), ]
    } else {
      data.list[[id]] <-
        obs[[id]][, c("x", "y", "time", "ID", "state", "group")]
    }
  }

  # flatten data
  data.mat <- do.call("rbind", data.list)
  data.mat <-
    data.mat[, c("x", "y", "time", "ID", "state", "group")]

  # initialise Hmat (rows of 0s for transitions)
  HmatAll <- matrix(0, nrow(data.mat), 4)
  HmatAll[which(!is.na(data.mat[, "x"])), ] <- Hmat

  # initial likelihood, logprior
  if (randomEffects) {
    oldllk <- do.call("sum", lapply(ids, function(id) {
      kalman_rcpp(
        data = data.mat[which(data.mat[, "ID"] == id), ],
        nbStates = nbStates,
        param = param[[id]],
        fixmu = mu,
        Hmat = HmatAll
        # normally distributed random effects
      )$llk + sum(dnorm(param[[id]][names(hyperparam_mu)][which(is.na(unlist(fixPar)))], hyperparam_mu[which(is.na(unlist(fixPar)))], hyperparam_sd[which(is.na(unlist(fixPar)))], log = TRUE))
    }))
    oldlogprior <- do.call("sum", lapply(ids, function(id) {
      getLogPrior(
        NULL,
        mu,
        fixPar,
        fixMu,
        priorFunc,
        priorArgs,
        rateparam,
        ratePriorFunc,
        ratePriorArgs
      )
    })) + sum(unlist(c(
      sapply(names(hyperparam_mu), FUN = function(par) {
        if (is.na(unlist(fixPar)[par])) {
          return(do.call(priorFunc[[par]], c(priorArgs[[par]], log = TRUE, x = hyperparam_mu[[par]])))
        }
      }),
      sapply(names(hyperparam_sd), FUN = function(par) {
        if (is.na(unlist(fixPar)[str_replace(par, "_sd", "")])) {
          return(do.call(priorFunc[[par]], c(priorArgs[[par]], log = TRUE, x = hyperparam_sd[[par]])))
        }
      })
    )))
  } else {
    oldllk <-
      kalman_rcpp(
        data = data.mat,
        nbStates = nbStates,
        param = param,
        fixmu = mu,
        Hmat = HmatAll
      )$llk
    oldlogprior <- getLogPrior(
      param,
      mu,
      fixPar,
      fixMu,
      priorFunc,
      priorArgs,
      rateparam,
      ratePriorFunc,
      ratePriorArgs
    )
  }

  ###############################
  ## Loop over MCMC iterations ##
  ###############################

  if (randomEffects) {
    allParam <- matrix(NA, nrow = nbIter - burnin, ncol = 2 * nbStates + 2 * nbHyperParam * nbStates)
    colnames(allParam) <- c(
      paste(
        c("mu_x[", "mu_y["),
        rep(1:nbStates, each = 2),
        c("]", "]"),
        sep = ""
      ),
      paste("tau_pos_mu[", 1:nbStates, "]", sep = ""),
      paste("tau_vel_mu[", 1:nbStates, "]", sep = ""),
      switch((nbParam == 7) + 1,
        paste("sigma_mu[", 1:nbStates, "]", sep = ""),
        paste(
          c("sigma_xx_mu[", "sigma_yy_mu[", "sigma_xy_mu["),
          rep(1:nbStates, each = 3),
          rep("]", nbStates * 3),
          sep = ""
        )
      ),
      paste("tau_pos_sd[", 1:nbStates, "]", sep = ""),
      paste("tau_vel_sd[", 1:nbStates, "]", sep = ""),
      switch((nbParam == 7) + 1,
        paste("sigma_sd[", 1:nbStates, "]", sep = ""),
        paste(
          c("sigma_xx_Sd[", "sigma_yy_sd[", "sigma_xy_sd["),
          rep(1:nbStates, each = 3),
          rep("]", nbStates * 3),
          sep = ""
        )
      )
    )
  } else {
    allParam <- matrix(NA, nrow = nbIter - burnin, ncol = nbParam * nbStates)
    colnames(allParam) <- c(
      paste("tau_pos[", 1:nbStates, "]", sep = ""),
      paste("tau_vel[", 1:nbStates, "]", sep = ""),
      switch((nbParam == 7) + 1,
        paste("sigma[", 1:nbStates, "]", sep = ""),
        paste(
          c("sigma_xx[", "sigma_yy[", "sigma_xy["),
          rep(1:nbStates, each = 3),
          rep("]", nbStates * 3),
          sep = ""
        )
      ),
      paste(
        c("mu_x[", "mu_y["),
        rep(1:nbStates, each = 2),
        c("]", "]"),
        sep = ""
      )
    )
  }
  row.names(allParam) <- format((burnin + 1):nbIter, scientific = FALSE, trim = TRUE)
  allLLk <- rep(NA, nbIter - burnin)
  names(allLLk) <- format((burnin + 1):nbIter, scientific = FALSE, trim = TRUE)
  # accParam <- rep(0, nbIter - burnin)
  accParam <- 0

  if (updateState) {
    idx <- seq(0, nbIter, thinStates)[which(seq(0, nbIter, thinStates) > burnin)]
    allStates <-
      matrix(NA, nrow = length(idx), ncol = nrow(track)) # uses a lot of memory if not thinning!
    row.names(allStates) <- format(idx, scientific = FALSE, trim = TRUE)
    rm(idx)
    # allLen <- matrix(NA, nrow = nbIter - burnin, ncol = length(ids))
    # accSwitch <- rep(0, nbIter - burnin)
    accSwitch <- 0
    if (!is.null(Q)) {
      allRates <-
        array(
          NA,
          dim = c(nbIter - burnin, nbStates * (nbStates - 1), length(ids)),
          dimnames = list(format((burnin + 1):nbIter, scientific = FALSE, trim = TRUE), NULL, NULL)
        )
    } else {
      allRateParam <- matrix(NA, nrow = nbIter - burnin, ncol = length(rateparam))
      row.names(allRateParam) <- format((burnin + 1):nbIter, scientific = FALSE, trim = TRUE)
    }
  }

  t0 <- Sys.time()
  for (iter in 1:nbIter) {
    if (iter < 100 || iter %% 100 == 0) {
      cat(
        "\33[2K\rIteration ",
        iter,
        "/",
        nbIter,
        "... ",
        vague_dt(
          difftime(Sys.time(), t0, units = "secs") / iter * (nbIter - iter),
          "short"
        ),
        " remaining (est)",
        ifelse(
          exists("accSwitch"),
          paste0(" -- accSwitch = ", round(accSwitch / iter * 100), "%"),
          ""
        ),
        " -- accParam = ",
        round(accParam / iter * 100),
        "%",
        sep = ""
      )
    }

    ################################################
    ## 1. Update discrete state process and rates ##
    ################################################

    if (updateState) {
      newData.list <- data.list
      newSwitch <- switch

      if (is.null(Q) && !is.na(model)) {
        # propose new rate params
        newRateParams <- proposeParams(rateparam, rateS)
      } else {
        newRateParams <- NULL
      }

      if (!is.null(newRateParams) && any(!constraintsFromDensityArgs(newRateParams[[2]], ratePriorArgs))) {
        acceptProb <- 0
      } else {
        acceptProb <- 0
        # catch autoreject "error"
        if (randomEffects) {
          upState <- lapply(ids, function(id) {
            try(updateState(
              obs = obs[[id]],
              nbStates = nbStates,
              knownStates = known[[id]],
              switch = switch[[id]
              ],
              updateLim = updateLim[[id]],
              param = param[[id]],
              mu = unlist(mu),
              Hmat = HmatAll[which(data.mat[, "ID"] == id), ],
              updateProbs = updateProbs[[id]],
              Q = switch(is.null(Q) + 1,
                Q[[id]],
                NULL
              ),
              # https://www.r-bloggers.com/2017/02/use-switch-instead-of-ifelse-to-return-a-null/
              rateparam = newRateParams[[2]],
              kappa = kappa,
              model = model
            ), !debug)
          })
        } else {
          upState <- lapply(ids, function(id) {
            try(updateState(
              obs = obs[[id]],
              nbStates = nbStates,
              knownStates = known[[id]],
              switch = switch[[id]
              ],
              updateLim = updateLim[[id]],
              param = param,
              mu = unlist(mu),
              Hmat = HmatAll[which(data.mat[, "ID"] == id), ],
              updateProbs = updateProbs[[id]],
              Q = switch(is.null(Q) + 1,
                Q[[id]],
                NULL
              ),
              # https://www.r-bloggers.com/2017/02/use-switch-instead-of-ifelse-to-return-a-null/
              rateparam = newRateParams[[2]],
              kappa = kappa,
              model = model
            ), !debug)
          })
        }
        if (all(sapply(upState, function(x) !inherits(x, "try-error")))) {
          newData.list <- lapply(ids, function(id) upState[[id]]$newData)
          newSwitch <- lapply(ids, function(id) upState[[id]]$newSwitch)

          # flatten data
          newData.mat <- do.call("rbind", newData.list)

          # update Hmat (rows of 0s for transitions)
          newHmatAll <- matrix(0, nrow(newData.mat), 4)
          newHmatAll[which(!is.na(newData.mat[, "x"])), ] <- Hmat

          # Calculate acceptance ratio
          if (randomEffects) {
            newllk <- do.call("sum", lapply(ids, function(id) {
              kalman_rcpp(
                data = newData.mat[which(newData.mat[, "ID"] == id), ],
                nbStates = nbStates,
                param = param[[id]],
                fixmu = unlist(mu),
                Hmat = newHmatAll
                # normally distributed random effects
              )$llk + sum(dnorm(param[[id]][names(hyperparam_mu)][which(is.na(unlist(fixPar)))], hyperparam_mu[which(is.na(unlist(fixPar)))], hyperparam_sd[which(is.na(unlist(fixPar)))], log = TRUE))
            }))
            newlogprior <- do.call("sum", lapply(ids, function(id) {
              getLogPrior(
                NULL,
                mu,
                fixPar,
                fixMu,
                priorFunc,
                priorArgs,
                newRateParams[[2]],
                ratePriorFunc,
                ratePriorArgs
              )
            })) + sum(unlist(c(
              sapply(names(hyperparam_mu), FUN = function(par) {
                if (is.na(unlist(fixPar)[par])) {
                  return(do.call(priorFunc[[par]], c(priorArgs[[par]], log = TRUE, x = hyperparam_mu[[par]])))
                }
              }),
              sapply(names(hyperparam_sd), FUN = function(par) {
                if (is.na(unlist(fixPar)[str_replace(par, "_sd", "")])) {
                  return(do.call(priorFunc[[par]], c(priorArgs[[par]], log = TRUE, x = hyperparam_sd[[par]])))
                }
              })
            )))
          } else {
            newllk <-
              kalman_rcpp(
                data = newData.mat,
                nbStates = nbStates,
                param = param,
                fixmu = unlist(mu),
                Hmat = newHmatAll
              )$llk
            newlogprior <- getLogPrior(
              param,
              mu,
              fixPar,
              fixMu,
              priorFunc,
              priorArgs,
              newRateParams[[2]],
              ratePriorFunc,
              ratePriorArgs
            )
          }
          logHR <- newllk + newlogprior - oldllk - oldlogprior
          acceptProb <-
            min(1, ifelse(is.na(exp(logHR)), 0, exp(logHR)))

          if (runif(1) < acceptProb) {
            # Accept new state sequence
            # accSwitch[iter] <- 1
            accSwitch <- accSwitch + 1
            switch <- newSwitch
            data.list <- newData.list
            obs <-
              lapply(data.list, function(data) {
                data[!is.na(data[, "x"]), ]
              })
            oldllk <- newllk
            rateparam <- newRateParams[[2]]
            oldlogprior <- newlogprior
            HmatAll <- newHmatAll
          }
        }
      }
      if (!is.na(model) && adapt && iter >= 100 && iter <= adapt) {
        newS <- adapt_S(rateS, newRateParams[[1]], acceptProb, iter)
        newS[is.na(newS)] <- rateS[is.na(newS)]
        rateS <- newS
      }

      if (!is.null(Q) && is.na(model)) {
        Q <- lapply(ids, function(id) {
          updateQ(
            nbStates = nbStates,
            data = data.list[[id]],
            switch = switch[[id]
            ],
            priorShape = priorShape,
            priorRate = priorRate,
            priorCon = priorCon,
            kappa = kappa
          )
        })
        names(Q) <- ids
      }
    }

    ###################################
    ## 2. Update movement parameters ##
    ###################################

    if (randomEffects) {
      newHyperParams <-
        proposeParams(c(hyperparam_mu, hyperparam_sd), S[seq_along(c(hyperparam_mu, hyperparam_sd)), seq_along(c(hyperparam_mu, hyperparam_sd))], nbStates)
      newParams <- lapply(ids, function(id) {
        newP <- proposeParams(param[[id]], randomS[[id]], nbStates)
        # overwrite fixed params
        newP[[2]][which(!is.na(unlist(fixPar)))] <- unlist(fixPar)[which(!is.na(unlist(fixPar)))]
        newP
      })
    } else {
      newParams <-
        proposeParams(param, S[seq_along(param), seq_along(param)], nbStates)
      # overwrite fixed params
      newParams[[2]][which(!is.na(unlist(fixPar)))] <- unlist(fixPar)[which(!is.na(unlist(fixPar)))]
    }

    if (randomEffects) {
      if (any(sapply(newParams, function(pars) {
        !constraintsFromDensityArgs(pars[[2]], priorArgs)
      })) || any(sapply(newParams, function(pars) {
        pars[[2]][1:nbStates] < pars[[2]][(nbStates + 1):(2 * nbStates)]
      })) || any(newHyperParams[[2]][which(str_detect(names(newHyperParams[[2]]), "_sd"))][which(is.na(unlist(fixPar)))] < 0)) {
        acceptProb <- 0
      } else {
        newMu <-
          proposeParams(mu, S[(2 * nbHyperParam * nbStates + 1):nrow(S), (2 * nbHyperParam * nbStates + 1):ncol(S)])
        # NB we could bound mu to -180,180 -90,90 with a different dist in proposeMus() but would need projected bounds
        newMu[[2]][which(!is.na(unlist(fixMu)))] <- unlist(fixMu)[which(!is.na(unlist(fixMu)))]
        # Calculate acceptance ratio
        data.mat <- do.call("rbind", data.list)
        data.mat <- data.mat[, c("x", "y", "time", "ID", "state")]
        newlogprior <- do.call("sum", lapply(ids, function(id) {
          getLogPrior(
            NULL,
            newMu[[2]],
            fixPar,
            fixMu,
            priorFunc,
            priorArgs,
            rateparam,
            ratePriorFunc,
            ratePriorArgs
          )
        })) + sum(unlist(c(
          sapply(names(newHyperParams[[2]]), FUN = function(par) {
            if (is.na(unlist(fixPar)[str_replace(par, "_sd", "")])) {
              return(do.call(priorFunc[[par]], c(priorArgs[[par]], log = TRUE, x = newHyperParams[[2]][[par]])))
            }
          })
        )))
        # TODO: could we partially accept proposals? this is absolute
        newllk <- do.call("sum", lapply(ids, function(id) {
          kalman_rcpp(
            data = data.mat[which(data.mat[, "ID"] == id), ],
            nbStates = nbStates,
            param = newParams[[id]][[2]],
            fixmu = newMu[[2]],
            Hmat = HmatAll
            # normally distributed random effects
            #TODO: subset to observed states?
          )$llk + sum(dnorm(newParams[[id]][[2]][names(hyperparam_mu)][which(is.na(unlist(fixPar)))], newHyperParams[[2]][which(!str_detect(names(newHyperParams[[2]]), "_sd"))][which(is.na(unlist(fixPar)))], newHyperParams[[2]][which(str_detect(names(newHyperParams[[2]]), "_sd"))][which(is.na(unlist(fixPar)))], log = TRUE))
        }))
        # mu <- as.vector(t(kalman$mu))
        logHR <- newllk + newlogprior - oldllk - oldlogprior
        acceptProb <- min(1, ifelse(is.na(exp(logHR)), 0, exp(logHR)))
      }

      if (runif(1) < acceptProb) {
        # Accept new parameter values
        # accParam[iter] <- 1
        accParam <- accParam + 1
        param <- lapply(ids, function(id) {
          newParams[[id]][[2]]
        })
        hyperparam_mu <- newHyperParams[[2]][which(!str_detect(names(newHyperParams[[2]]), "_sd"))]
        hyperparam_sd <- newHyperParams[[2]][which(str_detect(names(newHyperParams[[2]]), "_sd"))]
        mu <- newMu[[2]]
        oldllk <- newllk
        oldlogprior <- newlogprior
      }
    } else {
      if (any(!constraintsFromDensityArgs(newParams[[2]], priorArgs)) || any(newParams[[2]][1:nbStates] < newParams[[2]][(nbStates + 1):(2 * nbStates)])) {
        acceptProb <- 0
      } else {
        newMu <-
          proposeParams(mu, S[(length(param) + 1):nrow(S), (length(param) + 1):ncol(S)])
        # NB we could bound mu to -180,180 -90,90 with a different dist in proposeMus() but would need projected bounds
        newMu[[2]][which(!is.na(unlist(fixMu)))] <- unlist(fixMu)[which(!is.na(unlist(fixMu)))]

        # Calculate acceptance ratio
        data.mat <- do.call("rbind", data.list)
        data.mat <- data.mat[, c("x", "y", "time", "ID", "state")]
        newlogprior <- getLogPrior(
          newParams[[2]],
          newMu[[2]],
          fixPar,
          fixMu,
          priorFunc,
          priorArgs,
          rateparam,
          ratePriorFunc,
          ratePriorArgs
        )
        kalman <-
          kalman_rcpp(
            data = data.mat,
            nbStates = nbStates,
            param = newParams[[2]],
            fixmu = newMu[[2]],
            Hmat = HmatAll
          )
        newllk <- kalman$llk
        # mu <- as.vector(t(kalman$mu))
        logHR <- newllk + newlogprior - oldllk - oldlogprior
        acceptProb <- min(1, ifelse(is.na(exp(logHR)), 0, exp(logHR)))
      }

      if (runif(1) < acceptProb) {
        # Accept new parameter values
        # accParam[iter] <- 1
        accParam <- accParam + 1
        param <- newParams[[2]]
        mu <- newMu[[2]]
        oldllk <- newllk
        oldlogprior <- newlogprior
      }
    }

    # TODO: adapt and burnin redundant?
    if (adapt && iter >= 100 && iter <= adapt) {
      if (randomEffects) {
        if (any(is.na(unlist(fixPar)))) {
          newS <-
            adapt_S(S[seq_along(c(hyperparam_mu, hyperparam_sd)), seq_along(c(hyperparam_mu, hyperparam_sd))][is.na(c(unlist(fixPar), unlist(fixPar))), is.na(c(unlist(fixPar), unlist(fixPar)))], newHyperParams[[1]][is.na(c(unlist(fixPar), unlist(fixPar)))], acceptProb, iter)
          newS[is.na(newS)] <-
            S[seq_along(c(hyperparam_mu, hyperparam_sd)), seq_along(c(hyperparam_mu, hyperparam_sd))][is.na(c(unlist(fixPar), unlist(fixPar))), is.na(c(unlist(fixPar), unlist(fixPar)))][is.na(newS)]
          S[seq_along(c(hyperparam_mu, hyperparam_sd)), seq_along(c(hyperparam_mu, hyperparam_sd))][is.na(c(unlist(fixPar), unlist(fixPar))), is.na(c(unlist(fixPar), unlist(fixPar)))] <-
            newS
          for (id in ids) {
            newS <-
              adapt_S(randomS[[id]][seq_along(param[[id]]), seq_along(param[[id]])][is.na(unlist(fixPar)), is.na(unlist(fixPar))], newParams[[id]][[1]][is.na(unlist(fixPar))], acceptProb, iter)
            newS[is.na(newS)] <-
              randomS[[id]][seq_along(param[[id]]), seq_along(param[[id]])][is.na(unlist(fixPar)), is.na(unlist(fixPar))][is.na(newS)]
            randomS[[id]][seq_along(param[[id]]), seq_along(param[[id]])][is.na(unlist(fixPar)), is.na(unlist(fixPar))] <-
              newS
          }
        }

        muindex <-
          unique(cumsum(rep(
            !is.infinite(fixPar$tau_pos),
            each = 2
          )))

        if (any(is.na(unlist(fixMu)[muindex]))) {
          newS <-
            adapt_S(S[length(c(hyperparam_mu, hyperparam_sd)) + muindex, length(c(hyperparam_mu, hyperparam_sd)) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])], newMu[[1]][muindex][is.na(unlist(fixMu)[muindex])], acceptProb, iter)
          newS[is.na(newS)] <-
            S[length(c(hyperparam_mu, hyperparam_sd)) + muindex, length(c(hyperparam_mu, hyperparam_sd)) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])][is.na(newS)]
          S[length(c(hyperparam_mu, hyperparam_sd)) + muindex, length(c(hyperparam_mu, hyperparam_sd)) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])] <-
            newS
        }
      } else {
        if (any(is.na(unlist(fixPar)))) {
          newS <-
            adapt_S(S[seq_along(param), seq_along(param)][is.na(unlist(fixPar)), is.na(unlist(fixPar))], newParams[[1]][is.na(unlist(fixPar))], acceptProb, iter)
          newS[is.na(newS)] <-
            S[seq_along(param), seq_along(param)][is.na(unlist(fixPar)), is.na(unlist(fixPar))][is.na(newS)]
          S[seq_along(param), seq_along(param)][is.na(unlist(fixPar)), is.na(unlist(fixPar))] <-
            newS
        }

        muindex <-
          unique(cumsum(rep(
            !is.infinite(fixPar$tau_pos),
            each = 2
          )))

        if (any(is.na(unlist(fixMu)[muindex]))) {
          newS <-
            adapt_S(S[length(param) + muindex, length(param) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])], newMu[[1]][muindex][is.na(unlist(fixMu)[muindex])], acceptProb, iter)
          newS[is.na(newS)] <-
            S[length(param) + muindex, length(param) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])][is.na(newS)]
          S[length(param) + muindex, length(param) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])] <-
            newS
        }
      }
    }

    #########################
    ## Save posterior draw ##
    #########################
    if (iter > burnin) {
      if (randomEffects) {
        allParam[format(iter, scientific = FALSE, trim = TRUE), ] <-
          cbind(matrix(mu, ncol = 2 * nbStates), matrix(hyperparam_mu, ncol = length(hyperparam_mu)), matrix(hyperparam_sd, ncol = length(hyperparam_sd)))
      } else {
        allParam[format(iter, scientific = FALSE, trim = TRUE), ] <-
          cbind(matrix(param, ncol = length(param)), matrix(mu, ncol = 2 * nbStates))
      }
      if (updateState && !is.null(Q)) {
        allRates[format(iter, scientific = FALSE, trim = TRUE), , ] <-
          matrix(
            unlist(lapply(Q, function(q) {
              q[!diag(nbStates)]
            })),
            ncol = length(ids),
            nrow = nbStates * (nbStates - 1)
          )
      } else if (updateState) {
        allRateParam[format(iter, scientific = FALSE, trim = TRUE), ] <- rateparam
      }

      if (iter %% thinStates == 0 && updateState) {
        allStates[format(iter, scientific = FALSE, trim = TRUE), ] <-
          unlist(lapply(obs, function(ob) {
            ob[, "state"]
          }))
      }
      allLLk[format(iter, scientific = FALSE, trim = TRUE)] <- oldllk
    }
  }
  cat("\n")
  cat("Elapsed: ", pretty_dt(difftime(Sys.time(), t0, units = "secs")), sep = "")
  cat("\n")

  return(
    c(
      list(inits = inits),
      list(priors = priors),
      list(allParam = allParam),
      if (exists("allRates")) {
        list(allRates = allRates)
      },
      if (exists("allRateParam")) {
        list(allRateParam = allRateParam)
      },
      if (exists("allStates")) {
        list(allStates = allStates)
      },
      # if (exists("accSwitch")) {
      #  list(accSwitch = accSwitch)
      # },
      # list(accParam = accParam),
      # if (exists("allLen")) {
      #  list(allLen = allLen)
      # },
      list(allLLk = allLLk)
    )
  )
}

proposeParams <- function(param, S, nbStates) {
  u <- rnorm(length(param))
  thetasprime <- param + as.vector(S %*% u)
  return(list(u, thetasprime))
}

getLogPrior <- function(
    param, mu, fixPar, fixMu, priorFunc, priorArgs,
    rateparam, ratePriorFunc, ratePriorArgs) {
  return(
    sum(
      unlist(c(
        sapply(names(param), FUN = function(par) {
          if (is.na(unlist(fixPar)[par])) {
            return(do.call(priorFunc[[par]], c(priorArgs[[par]], log = TRUE, x = param[[par]])))
          }
        }),
        sapply(seq_len(length(mu)), FUN = function(i) {
          if (is.na(unlist(fixMu)[i]) & !is.na(mu[i])) {
            return(do.call(priorFunc[[paste0("mu", i)]], c(priorArgs[[paste0("mu", i)]], log = TRUE, x = mu[[i]])))
          }
        }),
        sapply(names(rateparam), FUN = function(par) {
          return(do.call(ratePriorFunc[[par]], c(ratePriorArgs[[par]], log = TRUE, x = rateparam[[par]])))
        })
      ))
    )
  )
}

constraintsFromDensityArgs <- function(param, args) {
  return(sapply(names(param), FUN = function(par) {
    if (is.na(param[par])) {
      return(TRUE)
    }
    lower <- ifelse(any(names(args[[par]]) %in% c("min", "lower")), param[par] >= args[[par]][[which(names(args[[par]]) %in% c("min", "lower"))]], TRUE)
    upper <- ifelse(any(names(args[[par]]) %in% c("max", "upper")), param[par] <= args[[par]][[which(names(args[[par]]) %in% c("max", "upper"))]], TRUE)
    return(lower && upper)
  }))
}

stopIfNotDataFrame <- function(df) {
  argname <- deparse(substitute(df))
  if (!is.data.frame(df)) {
    stop(paste0("argument '", argname, "' is not a data.frame"))
  }
}

isMissing <- function(obj, objs) {
  !sapply(objs, exists, where = obj)
}

stopIfMissingColumns <- function(df, cols) {
  argname <- deparse(substitute(df))
  missing <- df |> isMissing(cols)
  if (any(!(cols %in% colnames(df)))) {
    stop(
      paste0("argument '", argname, "' is missing required column(s): "),
      paste(cols[which(missing)], collapse = ", ")
    )
  }
}

stopIfMissingObjects <- function(l, objs) {
  argname <- deparse(substitute(l))
  missing <- l |> isMissing(objs)
  if (any(missing)) {
    stop(
      paste0("argument '", argname, "' is missing required object(s): "),
      paste(objs[which(missing)], collapse = ", ")
    )
  }
}

stopIfLengthNotEquals <- function(l, objs, lens) {
  argname <- deparse(substitute(l))
  for (obj in objs) {
    if (all(length(l[[obj]]) != lens)) {
      stop(
        paste0("argument '", argname, "$", obj, "'"),
        " has the wrong length, expected ",
        paste(lens, collapse = " or "),
        " but got ",
        length(l[[obj]])
      )
    }
  }
}

stopIfNestedLengthNotEquals <- function(l, obj, len) {
  argname <- deparse(substitute(l))
  if (all(sapply(l[[obj]], length) != len)) {
    stop(
      paste0("argument '", argname, "$", obj, "' has the wrong length, expected "),
      paste(len, collapse = " "),
      " but got ",
      paste(sapply(l[[obj]], length), collapse = " ")
    )
  }
}

stopIfDimNotEquals <- function(l, obj, dims) {
  argname <- deparse(substitute(l))
  if (all(dim(l[[obj]]) != dims)) {
    stop(
      paste0("argument '", argname, "$", obj, "' has the wrong dimensions, expected "),
      paste(dims, collapse = " "),
      " but got ",
      paste(dim(l[[obj]]), collapse = " ")
    )
  }
}
