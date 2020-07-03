
#' Run MCMC iterations
#'
#' @param track Dataframe of data, with columns "x", "y", "time", and "ID"
#' @param nbStates Number of states
#' @param nbIter Number of iterations
#' @param fixpar Vector of fixed parameter values (tau_pos, tau_vel, sigma), with NA for parameters to estimate
#' @param inits List of initial parameters
#' (tau_pos, tau_vel, sigma, Q, state)
#' @param priors List of parameters of prior distributions, with components:
#' \itemize{
#'   \item{"mean":} Vector of means for normal priors on movement parameters, of length 3*nbStates
#'   \item{"sd":} Vector of standard deviations for normal priors on movement parameters, of length
#' 3*nbStates
#'   \item{"shape":} Vector of shapes of gamma priors for the transition rates
#'   \item{"rate":} Vector of rates of gamma priors for the transition rates
#'   \item{"con":} Vector of concentrations of Dirichlet priors for transition probabilities
#' }
#' @param props List of parameters of proposal distributions, with components:
#' \itemize{
#'   \item{"tau_posSD":} Scalar standard deviation for normal proposal distribution of tau_pos
#'   \item{"tau_velSD":} Scalar standard deviation for normal proposal distribution of tau_vel
#'   \item{"sigmaSD":} Scalar standard deviation for normal proposal distribution of sigma
#'   \item{"updateLim":} Vector of two values: min and max length of updated state sequence
#'   \item{"updateProbs":} Probability for each element of updateLim[1]:updateLim[2] (if NULL,
#'   all values are equiprobable)
#' }
#' @param tunes List of tuning parameters, with components:
#' \itemize{
#'   \item{"thinStates":} Thinning factor for the posterior state sequences (needed because
#' of memory limitations)
#' }
#' @param Hmat Matrix of observation error variance (four columns, and one row
#' for each row of data)
#'
#' @param updateState Logical. If FALSE, the state process is not updated
#' (for exploratory analysis only, useful for testing single state models)
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats dnorm runif rnorm rexp rgamma
#' @importFrom parallel mclapply
#' @importFrom prettyunits vague_dt pretty_dt
#' @export
#'
#' @useDynLib MSctmm
runMCMC <- function(track, nbStates, nbIter, fixpar=NULL, inits, priors, props, tunes, Hmat,
                    updateState=TRUE, mc.cores = 1)
{

    ######################
    ## Unpack arguments ##
    ######################
    # initial parameters
    tau_pos <- inits$tau_pos
    tau_vel <- inits$tau_vel
    #NB considering only isotropic case for now, this will eventually need to be a matrix
    #kalman_rcpp and MakeQ will need to accept a matrix
    sigma <- inits$sigma
    param <- c(tau_pos, tau_vel, sigma)
    Q <- rep(list(inits$Q), length( unique( track$ID ) ) )
    names(Q) <- unique( track$ID )
    state0 <- inits$state

    if( is.null( fixpar ) ) {
      fixpar <- list( tau_pos = rep(NA, nbStates),
                      tau_vel = rep(NA, nbStates),
                      sigma = rep(NA, nbStates))
    } else {
      fixpar <- fixpar
    }
    if(length(unlist(inits[1:3]))!=length(unlist(fixpar)))
      stop("'fixpar' has the wrong length")

    if(is.null(tau_pos) | is.null(tau_vel) | is.null(sigma) | is.null(Q) | is.null(state0) )
        stop("'inits' should have components tau_pos, tau_vel, sigma, Q, and state")

    # Kalman filter parameters
    if(is.null(Hmat))
      stop("'Hmat' should not be null")

    # unpack prior parameters
    priorMean <- priors$mean
    priorSD <- priors$sd
    priorShape <- priors$shape
    priorRate <- priors$rate
    priorCon <- priors$con
    if(is.null(priorMean) | is.null(priorSD) | is.null(priorShape) |
       is.null(priorRate) | is.null(priorCon)) {
        stop("'priors' should have components mean, sd, shape, rate, and con")
    }

    # unpack proposal parameters
    tau_posPropSD <- props$tau_posSD
    tau_velPropSD <- props$tau_velSD
    sigmaPropSD <- props$sigmaSD
    updateLim <- props$updateLim
    updateProbs <- props$updateProbs
    if(is.null(tau_posPropSD) | is.null(tau_velPropSD) | is.null(sigmaPropSD) | is.null(updateLim))
        stop("'props' should have components tau_posSD, tau_velSD, sigmaSD, and updateLim")

    if(is.null(updateProbs))
        updateProbs <- rep(1, length(updateLim[1]:updateLim[2]))/length(updateLim[1]:updateLim[2])

    if(length(updateLim[1]:updateLim[2])!=length(updateProbs))
        stop("'updateProbs' has the wrong length")

    # unpack tuning parameters
    thinStates <- tunes$thinStates

    ####################
    ## Initialisation ##
    ####################
    # Prepare data structures
    if( !( "ID" %in% colnames( track ) ) ) {
      warning( "'track' should have column 'ID', assuming data is one individual"  )
      track$ID <- 1
    }

    ids <- as.character(unique(track$ID))
    obs <- list()
    switch <- list()
    data <- list()
    HmatAll <- list()

    for( i in ids ) {
      nbObs <- nrow(track[ which( track$ID == i ), ])
      obs[[ i ]] <- matrix( c( track[ which( track$ID == i ), "x" ],
                               track[ which( track$ID == i ), "y" ],
                               track[ which( track$ID == i ), "time" ],
                               state0[ which( track$ID == i )  ] ),
                            ncol = 4 )
      colnames( obs[[ i ]] ) <- c( "x", "y", "time", "state" )
      indSwitch <- which( obs[[ i ]][ -1, "state" ] != obs[[ i ]][ -nbObs, "state" ] ) + 1
      switch[[ i ]] <- matrix( c( obs[[ i ]][ indSwitch, "time" ] - 0.001, rle( obs[[ i ]][ , "state" ] )$values[ -1 ] ), ncol = 2 )
      colnames( switch[[ i ]] ) <- c( "time", "state" )
      if( nrow( switch[[ i ]] ) ) {
        data[[ i ]] <- rbind( obs[[ i ]] ,cbind( NA, NA, switch[[ i ]] ) )
      } else {
        data [[ i ]] <- obs[[ i ]]
      }
      data[[ i ]] <- data[[ i ]][ order( data[[ i ]][ , "time" ] ), ]

      # initialise Hmat (rows of 0s for transitions)
      HmatAll[[ i ]] <- matrix( 0, nrow( data[[ i ]] ), 4 )
      HmatAll[[ i ]][ which( !is.na( data[[ i ]][ , "x" ] ) ), ] <- Hmat[ which( track$ID == i ), ]
    }

    # initial likelihood
    oldllk <- mclapply( ids, function( id ) { kalman_rcpp( data = data[[ id ]], param = param, Hmat = HmatAll[[ id ]] ) }, mc.cores = mc.cores )
    names(oldllk) <- ids

    # initial log-prior
    oldlogprior <- sum( dnorm( log( param[ unlist(is.na(fixpar)) ] ), priorMean[ unlist(is.na(fixpar)) ], priorSD[ unlist(is.na(fixpar)) ], log = TRUE ) )

    ###############################
    ## Loop over MCMC iterations ##
    ###############################
    allparam <- matrix(NA,nrow=nbIter,ncol=3*nbStates)
    allrates <- array( NA, dim = c( nbIter, nbStates*(nbStates-1), length( ids ) ) )
    allstates <- matrix( NA, nrow = nbIter / thinStates, ncol = nrow( track ) ) # uses a lot of memory!
    accSwitch <- rep(0,nbIter)
    accParam <- rep(0,nbIter)
    allLen <- rep(NA,nbIter)

    t0 <- Sys.time()
    for(iter in 1:nbIter) {
        if(iter%%100==0)
            cat("\rIteration ", iter, "/", nbIter, "... ",
                vague_dt(difftime(Sys.time(),t0,units="secs")/iter*(nbIter-iter),"short"),
                " remaining (est)",
                " -- accSwitch = ", round(sum(accSwitch)/iter*100), "%",
                " -- accParam = ", round(sum(accParam)/iter*100), "%", sep="")

        ######################################
        ## 1. Update discrete state process ##
        ######################################
        if(updateState) {

            # pick an individual at random for which we will update the state sequence
            id <- sample( ids, 1 )
            upState <- updateState( obs = obs[[ id ]], switch = switch[[ id ]], updateLim = updateLim,
                                   updateProbs = updateProbs, Q = Q[[ id ]])
            newData <- data
            newData[[ id ]] <- upState$newData
            newSwitch <- switch
            newSwitch[[ id ]] <- upState$newSwitch
            allLen[iter] <- upState$len

            # update Hmat (rows of 0s for transitions)
            newHmatAll <- HmatAll
            newHmatAll[[ id ]] <- matrix( 0, nrow( newData[[ id ]] ), 4 )
            newHmatAll[[ id ]][ which( !is.na( newData[[ id ]] [ , "x" ] ) ), ] <- Hmat[ which( track$ID == id ), ]

            # Calculate acceptance ratio
            newllk <- oldllk
            newllk[[ id ]] <- kalman_rcpp( data = newData[[ id ]], param = param, Hmat = newHmatAll[[ id ]] )
            logHR <- do.call( 'sum', newllk ) - do.call( 'sum', oldllk )

            if(log(runif(1))<logHR) {
                # Accept new state sequence
                accSwitch[iter] <- 1
                switch <- newSwitch
                data <- newData
                obs <- lapply(data,function(data) {data[!is.na(data[,"x"]),]})
                oldllk <- newllk
                HmatAll <- newHmatAll
            }
        }

        ###################################
        ## 2. Update movement parameters ##
        ###################################
        # On working scale
        tau_posprimeW <- rnorm(nbStates,log(tau_pos),tau_posPropSD)
        tau_velprimeW <- rnorm(nbStates,log(tau_vel),tau_velPropSD)
        sigmaprimeW <- rnorm(nbStates,log(sigma),sigmaPropSD)
        newlogprior <- sum(dnorm(c(tau_posprimeW,tau_velprimeW,sigmaprimeW)[ unlist(is.na(fixpar)) ],priorMean[ unlist(is.na(fixpar)) ],priorSD[ unlist(is.na(fixpar)) ],log=TRUE))

        # On natural scale
        tau_posprime <- vector( 'double',nbStates )
        tau_velprime <- vector( 'double',nbStates )
        sigmaprime <- vector( 'double',nbStates )

        tau_posprime[ is.na(fixpar$tau_pos) ] <- exp(tau_posprimeW)[ is.na(fixpar$tau_pos) ]
        tau_posprime[ !is.na(fixpar$tau_pos) ] <- fixpar$tau_pos[ !is.na(fixpar$tau_pos) ]
        tau_velprime[ is.na(fixpar$tau_vel) ] <- exp(tau_velprimeW)[ is.na(fixpar$tau_vel) ]
        tau_velprime[ !is.na(fixpar$tau_vel) ] <- fixpar$tau_vel[ !is.na(fixpar$tau_vel) ]
        sigmaprime[ is.na(fixpar$sigma) ] <- exp(sigmaprimeW)[ is.na(fixpar$sigma) ]
        sigmaprime[ !is.na(fixpar$sigma) ] <- fixpar$sigma[ !is.na(fixpar$sigma) ]

        # Calculate acceptance ratio
        newllk <- mclapply( ids, function( id ) { kalman_rcpp( data = data[[ id ]], param = c( tau_posprime, tau_velprime, sigmaprime ), Hmat = HmatAll[[ id ]] ) }, mc.cores = mc.cores )
        names(newllk) <- ids
        logHR <- do.call( 'sum', newllk ) + newlogprior - do.call( 'sum', oldllk ) - oldlogprior

        if(log(runif(1))<logHR) {
            # Accept new parameter values
            accParam[iter] <- 1
            tau_pos <- tau_posprime
            tau_vel <- tau_velprime
            sigma <- sigmaprime
            param <- c( tau_pos, tau_vel, sigma )
            oldllk <- newllk
            oldlogprior <- newlogprior
        }

        ###############################
        ## 3. Update switching rates ##
        ###############################
        Q <- lapply( ids, function( id ) { updateQ( nbStates = nbStates, data = data[[ id ]], switch = switch[[ id ]],
                     priorShape = priorShape, priorRate = priorRate,
                     priorCon = priorCon ) } )
        names(Q) <- ids

        #########################
        ## Save posterior draw ##
        #########################
        allparam[iter,] <- param
        allrates[ iter, , ] <- matrix( unlist( lapply( Q, function( q ){ q[ !diag( nbStates ) ] } ) ), ncol = length( ids ), nrow = nbStates*(nbStates-1) )
        if(iter%%thinStates==0)
            allstates[iter/thinStates,] <- unlist( lapply( obs, function( ob ) { ob[ , "state" ] } ) )
    }
    cat("\n")
    cat("Elapsed: ", pretty_dt(difftime(Sys.time(),t0,units="secs")),sep="")
    cat("\n")

    return(list(allparam = allparam,
                allrates = allrates,
                allstates = allstates,
                accSwitch = accSwitch,
                allLen = allLen))
}
