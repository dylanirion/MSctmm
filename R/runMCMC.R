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
#'   * `state`: vector. Initial state sequence, length `nrow(track)`
#'   * `rateparam`: vector. Length dependent on model choice
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
#'   * `func`: list of distribution function names, of length `5 * nbStates` or `7 * nbStates` for movement parameters + `length(rateparam)` when estimating rate parameters `rateparam` by MH
#'   * `args`: list of list containg function args for distribution functions, of length `5 * nbStates` or `7 * nbStates` for movement parameters + `length(rateparam)` when estimating rate parameters `rateparam` by MH
#'   * `shape`: vector. Shapes of gamma priors for the transition rates when Gibbs sampling
#'   * `rate`: vector. Rates of gamma priors for the transition rates when Gibbs sampling
#'   * `con`: vector. Concentrations of Dirichlet priors for transition probabilities when Gibbs sampling
#' @param props list. Parameters of proposal distributions, with components:
#'   * `S`: matrix. Initial value for the lower triangular matrix of RAM algorithm, so that the covariance matrix of the proposal distribution is `SS'`.
#'   Dimensions `(5 * nbStates, 5 * nbStates)` when not modelling rate
#'   parameters (`model` is `NA`) and `(5 * nbStates + length(rateparam), 5 * nbStates + length(rateparam))`
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
runMCMC <- function(track,
                    nbStates,
                    nbIter,
                    burnin,
                    inits,
                    fixed,
                    priors,
                    props,
                    tunes,
                    Hmat,
                    updateState = TRUE,
                    adapt = FALSE,
                    model = NA,
                    debug = FALSE,
                    silent = FALSE) {
  #TODO: a chainable MCMC sampler builder would be cool
  # model = new Model()
  # mod.addStep()
  options(warn = 1) # show warnings as they occur
  # Check track df
  if (!is.data.frame(track)) {
    stop("argument 'track' is not a data.frame")
  }
  if (any(!(c("x", "y", "time") %in% colnames(track)))) {
    stop(
      "argument 'track' is missing required column(s): ",
      paste(c("x", "y", "time")[which(!(c("x", "y", "time") %in% colnames(track)))], collapse = ", ")
    )
  }

  if (!("ID" %in% colnames(track))) {
    warning("argument 'track' should have column 'ID', assuming data is one individual")
    track$ID <- 1
  } else {
    track$ID <- as.numeric(as.factor(track$ID))
  }

  # Check inits arguments and lengths
  if (is.null(inits$tau_pos) ||
    is.null(inits$tau_vel) ||
    is.null(inits$sigma) || is.null(inits$mu) || is.null(inits$state)) {
    stop("argument 'inits' missing: ", paste(c(
      "tau_pos", "tau_vel", "sigma", "mu", "state"
    )[which(sapply(inits[c("tau_pos", "tau_vel", "sigma", "mu", "state")], is.null))], collapse = ", "))
  }
  for (arg in c("tau_pos", "tau_vel")) {
    if (length(inits[[arg]]) != nbStates) {
      stop(
        "argument 'inits$",
        arg,
        "' has the wrong length, expected ",
        nbStates,
        " but got ",
        length(inits[[arg]])
      )
    }
  }
  if (length(inits$sigma) != nbStates &&
    length(inits$sigma) != 3 * nbStates) {
    stop(
      "argument 'inits$sigma has the wrong length, expected ",
      nbStates,
      " or ",
      3 * nbStates,
      " but got ",
      length(inits$sigma)
    )
  }
  if (all(sapply(inits$mu, length) != rep(2, nbStates))) {
    stop(
      "argument 'inits$mu' has the wrong length, expected ",
      paste(rep(2, nbStates), collapse = " "),
      " but got ",
      paste(sapply(inits$mu, length), collapse = " ")
    )
  }
  if (length(inits$state) != nrow(track)) {
    stop(
      "'inits$state' has the wrong length, expected ",
      nrow(track),
      " but got ",
      length(inits$state)
    )
  }
  if (updateState &&
    (
      is.null(fixed$kappa) ||
        (is.null(inits$rateparam) && !is.na(model))
    )) {
    stop(
      "updateState is TRUE, expected ",
      paste(
        c("fixed$kappa")[which(is.null(fixed$kappa))],
        c("inits$rateparam")[which(sapply(inits[c("rateparam")], \(param) is.null(param) && !is.na(model)))],
        collapse = ", "
      ),
      " to be specified"
    )
  }

  # Check fixed arguments and lengths
  for (arg in c("tau_pos", "tau_vel", "sigma")) {
    if (is.null(fixed[[arg]])) {
      fixed[[arg]] <- rep(NA, nbStates)
    }
    if (length(fixed[[arg]]) != length(inits[[arg]])) {
      stop(
        "argument 'fixed$",
        arg,
        "' has the wrong length, expected ",
        length(inits[[arg]]),
        " but got ",
        length(fixed[[arg]])
      )
    }
  }
  if (is.null(fixed$mu)) {
    fixed$mu <- rep(list(c(NA, NA)), nbStates)
  } else {
    if (all(sapply(fixed$mu, length) != rep(2, nbStates))) {
      stop(
        "argument 'fixed$mu' has the wrong length, expected ",
        paste(rep(2, nbStates), collapse = " "),
        " but got ",
        paste(sapply(fixed$mu, length), collapse = " ")
      )
    }
  }
  if (is.null(fixed$knownStates) || is.na(fixed$knownStates)) {
    fixed$knownStates <- rep(NA, nrow(track))
  } else {
    if (length(fixed$knownStates) != nrow(track)) {
      stop(
        "argument 'fixed$knownStates' has the wrong length, expected ",
        nrow(track),
        " but got ",
        length(fixed$knownStates)
      )
    }
  }
  if (length(fixed$kappa) != 1) {
    stop(
      "argument 'fixed$kappa' has the wrong length, expected ",
      1,
      " but got ",
      length(fixed$kappa)
    )
  }

  # movement parameters (including mu)
  nbParam <- 4 + length(inits$sigma) / nbStates

  # Check priors arguments and lengths
  # TODO smarter calculation when args are fixed
  if (is.null(priors$func) || is.null(priors$args)) {
    stop("argument 'priors' missing: ", paste(c("func", "args")[which(sapply(priors[c("func", "args")], is.null))], collapse = ", "))
  }
  for (arg in c("func", "args")) {
    if (is.na(model) && length(priors[[arg]]) != nbParam * nbStates) {
      stop(
        "argument 'priors$",
        arg,
        "' has the wrong length, expected ",
        nbParam * nbStates,
        " but got ",
        length(priors[[arg]])
      )
    }
    if (!is.na(model) &&
      length(priors[[arg]]) != nbParam * nbStates + length(inits$rateparam)) {
      stop(
        "argument 'priors$",
        arg,
        "' has the wrong length, expected ",
        nbParam * nbStates + length(inits$rateparam),
        " but got ",
        length(priors[[arg]])
      )
    }
  }
  if (is.na(model) &&
    (is.null(priors$shape) ||
      is.null(priors$rate) || is.null(priors$con))) {
    stop(
      "argument 'model' is NA, but argument 'priors' missing: ",
      paste(c("shape", "rate", "con")[which(sapply(priors[c("shape", "rate", "con")], is.null))], collapse = ", ")
    )
  }
  if (is.na(model)) {
    for (arg in c("shape", "rate", "con")) {
      if (length(priors[[arg]]) != nbStates) {
        stop(
          "argument 'priors$",
          arg,
          "' has the wrong length, expected ",
          nbStates,
          " but got ",
          length(priors[[arg]])
        )
      }
    }
  }

  # Check props arguments and lengths
  if (is.null(props$S) || (updateState && is.null(props$updateLim))) {
    stop("argument 'props' missing: ", paste(c("S")[which(is.null(props$S))], c("S")[which(updateState &
      is.null(props$S))], collapse = ", "))
  }
  if (is.na(model) &&
    all(dim(props$S) != c(nbParam * nbStates, nbParam * nbStates))) {
    stop(
      "argument 'props$S' has the wrong dimensions, expected ",
      paste(c(nbParam * nbStates, nbParam * nbStates), collapse = ", "),
      " but got ",
      paste(dim(props$S), collapse = ", ")
    )
  }
  if (!is.na(model) &&
    all(
      dim(props$S) != c(
        nbParam * nbStates + length(inits$rateparam),
        nbParam * nbStates + length(inits$rateparam)
      )
    )) {
    stop(
      "argument 'props$S' has the wrong dimensions, expected ",
      paste(
        c(
          nbParam * nbStates + length(inits$rateparam),
          nbParam * nbStates + length(inits$rateparam)
        ),
        collpase = ", "
      ),
      " but got ",
      paste(dim(props$S), collapse = ", ")
    )
  }
  if (updateState && "list" %in% class(props$updateLim)) {
    if (!all(sapply(props$updateLim, length) == rep(2, length(unique(track$ID))))) {
      stop(
        "argument 'props$updateLim' has the wrong length, expected ",
        paste(rep(2, length(
          unique(track$ID)
        )), collapse = " "),
        " but got ",
        paste(sapply(props$updateLim, length), collapse = " ")
      )
    }
  }
  if (updateState &&
    !"list" %in% class(props$updateLim) && length(props$updateLim) != 2) {
    stop(
      "argument 'props$updateLim' has the wrong length, expected ",
      2,
      " but got ",
      length(props$updateLim)
    )
  }
  if (updateState &&
    !is.null(props$updateProbs) &&
    class(props$updateLim) != class(props$updateProbs)) {
    stop(
      "argument 'props$updateProbs' provided but of wrong type, expected ",
      class(props$updateLim),
      " but got ",
      class(props$updateProbs)
    )
  }
  if (updateState &&
    !is.null(props$updateProbs) &&
    "list" %in% class(props$updateProbs) &&
    all(sapply(props$updateLim, length) != sapply(props$updateProbs, length))) {
    stop(
      "argument 'props$updateProbs' has the wrong length, expected ",
      paste(sapply(props$updateLim, length), collapse = " "),
      " but got ",
      paste(sapply(props$updateProbs, length), collapse = " ")
    )
  }
  if (updateState &&
    !is.null(props$updateProbs) &&
    !"list" %in% class(props$updateProbs) &&
    length(props$updateProbs) != 2) {
    stop(
      "argument 'props$updateProbs' has the wrong length, expected ",
      2,
      " but got ",
      length(props$updateProbs)
    )
  }

  if (!updateState && !is.na(model) && ((!is.null(priors$shape) || !is.null(priors$rate)))) {
    warning(
      "argument 'updateState' is FALSE but 'model' is not NA, ignoring ",
      paste(c(
        "priors$shape", "priors$rate"
      )[which(!sapply(priors[c("shape", "rate")], is.null))], collapse = ", ")
    )
  }

  if (is.na(model) &&
    all(dim(props$S) > c(nbParam * nbStates, nbParam * nbStates))) {
    warning("argument 'model' is NA, ignoring extra dimensions of 'props$S")
  }

  # TODO check any NA in Hmat

  ######################
  ## Unpack arguments ##
  ######################

  # unpack fixed parameters
  fixPar <- fixed[c("tau_pos", "tau_vel", "sigma")]
  fixMu <- fixed$mu
  knownStates <- fixed$knownStates

  # initial movement parameters
  param <- unlist(inits[c("tau_pos", "tau_vel", "sigma")])
  mu <- unlist(inits$mu)
  state <- inits$state

  # overwrite fixed params
  param[which(!is.na(unlist(fixPar)))] <-
    unlist(fixPar)[which(!is.na(unlist(fixPar)))]
  mu[which(!is.na(unlist(fixMu)))] <-
    unlist(fixMu)[which(!is.na(unlist(fixMu)))]
  state[which(!is.na(knownStates))] <-
    knownStates[which(!is.na(knownStates))]

  # ensure tau_p > tau_v
  param[1:(2 * nbStates)] <- c(
    pmax(param[1:nbStates], param[(nbStates + 1):(2 * nbStates)]),
    pmin(param[1:nbStates], param[(nbStates + 1):(2 * nbStates)])
  )

  # unpack movement priors
  priorFunc <- priors$func[1:(nbParam * nbStates)]
  priorArgs <- priors$args[1:(nbParam * nbStates)]

  # unpack rate priors
  if (!updateState || is.na(model)) {
    ratePriorFunc <- NULL
    ratePriorArgs <- NULL
    rateparam <- NULL
    priorShape <- priors$shape
    priorRate <- priors$rate
  } else {
    rateS <-
      props$S[(nbParam * nbStates + 1):nrow(props$S), (nbParam * nbStates + 1):ncol(props$S)]
    ratePriorFunc <-
      priors$func[(nbParam * nbStates + 1):length(priors$func)]
    names(ratePriorFunc) <- sub("^rateparam\\.", "", names(ratePriorFunc))
    ratePriorArgs <-
      priors$args[(nbParam * nbStates + 1):length(priors$args)]
    names(ratePriorArgs) <- sub("^rateparam\\.", "", names(ratePriorArgs))
    rateparam <- inits$rateparam
    priorShape <- NULL
    priorRate <- NULL
  }
  priorCon <- priors$con
  kappa <- fixed$kappa

  # unpack proposal parameters
  S <- props$S[1:(nbParam * nbStates), 1:(nbParam * nbStates)]
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
    switch[[id]] <- cbind(
      "time" = obs[[id]][indSwitch, "time"] - 0.001,
      "state" = rle(obs[[id]][, "state"])$values[-1]
    )
    # colnames(switch[[id]]) <- c("time", "state")

    if (!all(is.na(switch[[id]]))) {
      data.list[[id]] <- rbind(
        obs[[id]][, c("x", "y", "time", "ID", "state", "group")],
        cbind(
          "x" = NA,
          "y" = NA,
          "time" = switch[[id]][, "time"],
          "ID" = id,
          "state" = switch[[id]][, "state"],
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
  # initial likelihood
  oldllk <-
    kalman_rcpp(
      data = data.mat,
      nbStates = nbStates,
      param = param,
      fixmu = mu,
      Hmat = HmatAll
    )$llk

  # initial log-prior
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

  ###############################
  ## Loop over MCMC iterations ##
  ###############################

  allParam <- matrix(NA, nrow = nbIter - burnin, ncol = nbParam * nbStates)
  row.names(allParam) <- format((burnin + 1):nbIter, scientific = FALSE, trim = TRUE)
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
    if (!is.na(model)) {
      allRateParam <- matrix(NA, nrow = nbIter - burnin, ncol = length(rateparam))
      row.names(allRateParam) <- format((burnin + 1):nbIter, scientific = FALSE, trim = TRUE)
      colnames(allRateParam) <- names(inits$rateparam)
    }
  }

  t0 <- Sys.time()
  for (iter in 1:nbIter) {
    if ((iter < 100 || iter %% 100 == 0) && !silent) {
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

      if (!is.na(model)) {
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
        upState <- lapply(ids, function(id) {
          try(updateState(
            obs = obs[[id]],
            nbStates = nbStates,
            knownStates = known[[id]],
            switch = switch[[id]],
            updateLim = updateLim[[id]],
            param = param,
            mu = unlist(mu),
            Hmat = HmatAll[which(data.mat[, "ID"] == id), ],
            updateProbs = updateProbs[[id]],
            rateparam = if (is.na(model)) {
              list(
                data = rbind(
                  data.list[[id]][1, c("time", "state")], # to include first interval
                  switch[[id]],
                  data.list[[id]][nrow(data.list[[id]]), c("time", "state")] # to include last interval
                ),
                priorShape = priorShape, priorRate = priorRate, priorCon = priorCon
              )
            } else {
              list(
                data = rbind(
                  data.list[[id]][1, c("time", "state")], # to include first interval
                  switch[[id]],
                  data.list[[id]][nrow(data.list[[id]]), c("time", "state")] # to include last interval
                ),
                params = newRateParams[[2]],
                priorCon = priorCon)
            },
            kappa = kappa,
            model = model
          ), !debug)
        })
        if (all(sapply(upState, function(x) !inherits(x, "try-error")))) {
          newData.list <- lapply(ids, function(id) upState[[id]]$newData)
          newSwitch <- lapply(ids, function(id) upState[[id]]$newSwitch)

          # flatten data
          newData.mat <- do.call("rbind", newData.list)

          # update Hmat (rows of 0s for transitions)
          newHmatAll <- matrix(0, nrow(newData.mat), 4)
          newHmatAll[which(!is.na(newData.mat[, "x"])), ] <- Hmat

          # Calculate acceptance ratio
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
      #if (!is.na(model) && adapt && iter >= 1000 && iter <= adapt) {
      #  newS <- adapt_S(rateS, newRateParams[[1]], acceptProb, iter)
      #  newS[is.na(newS)] <- rateS[is.na(newS)]
      #  rateS <- newS
      #}
    }

    ###################################
    ## 2. Update movement parameters ##
    ###################################

    newParams <-
      proposeParams(param, S[seq_along(param), seq_along(param)], nbStates)

    # overwrite fixed params
    newParams[[2]][which(!is.na(unlist(fixPar)))] <- unlist(fixPar)[which(!is.na(unlist(fixPar)))]

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

    # TODO: adapt and burnin redundant?
    if (adapt && iter >= 1000 && iter <= adapt) {
      if (any(is.na(unlist(fixPar)))) {
        newS <- adapt_S(S[seq_along(param), seq_along(param)][is.na(unlist(fixPar)), is.na(unlist(fixPar))], newParams[[1]][is.na(unlist(fixPar))], acceptProb, iter)
        newS[is.na(newS)] <- S[seq_along(param), seq_along(param)][is.na(unlist(fixPar)), is.na(unlist(fixPar))][is.na(newS)]
        S[seq_along(param), seq_along(param)][is.na(unlist(fixPar)), is.na(unlist(fixPar))] <- newS
      }


      muindex <-
        which(rep(
          !is.infinite(fixPar$tau_pos),
          each = 2
        ))
      if (any(is.na(unlist(fixMu)[muindex]))) {
        newS <- adapt_S(S[length(param) + muindex, length(param) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])], newMu[[1]][muindex][is.na(unlist(fixMu)[muindex])], acceptProb, iter)
        newS[is.na(newS)] <- S[length(param) + muindex, length(param) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])][is.na(newS)]
        S[length(param) + muindex, length(param) + muindex][is.na(unlist(fixMu)[muindex]), is.na(unlist(fixMu)[muindex])] <- newS
      }
    }

    #########################
    ## Save posterior draw ##
    #########################
    if (iter > burnin) {
      allParam[format(iter, scientific = FALSE, trim = TRUE), ] <-
        cbind(matrix(param, ncol = length(param)), matrix(mu, ncol = 2 * nbStates))
      if (updateState && !is.na(model)) {
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
  if (!silent) {
    cat("\n")
    cat("Elapsed: ", pretty_dt(difftime(Sys.time(), t0, units = "secs")), sep = "")
    cat("\n")

    if (debug) {
      cat(diag(S), "\n")
    }
  }

  return(
    c(
      list(inits = inits),
      list(priors = priors),
      list(allParam = allParam),
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

getLogPrior <- function(param,
                        mu,
                        fixPar,
                        fixMu,
                        priorFunc,
                        priorArgs,
                        rateparam,
                        ratePriorFunc,
                        ratePriorArgs) {
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
            return(do.call(priorFunc[[length(param) + i]], c(priorArgs[[length(param) + i]], log = TRUE, x = mu[[i]])))
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
