
#' Run MCMC iterations
#'
#' @param track Dataframe of data, with columns \code{"x"}, \code{"y"}, \code{"time"}, and \code{"ID"}
#' @param nbStates Number of states
#' @param nbIter Number of iterations
#' @param fixpar List of fixed parameter values (\code{"tau_pos"}, \code{"tau_vel"}, \code{"sigma"}), with \code{"NA"} for parameters to estimate
#' @param inits List of initial parameters (tau_pos, tau_vel, sigma, Q, state)
#' @param priors List of parameters of prior distributions, with components:
#' \itemize{
#'   \item{"mean":} Vector of means for normal priors on movement parameters, of length \code{"3*nbStates"}
#'   \item{"sd":} Vector of standard deviations for normal priors on movement parameters, of length \code{"3*nbStates"}
#'   \item{"shape":} Vector of shapes of gamma priors for the transition rates
#'   \item{"rate":} Vector of rates of gamma priors for the transition rates
#'   \item{"con":} Vector of concentrations of Dirichlet priors for transition probabilities
#' }
#' @param props List of parameters of proposal distributions, with components:
#' \itemize{
#'   \item{"S":} Initial value for the lower triangular matrix of RAM algorithm, so that the covariance matrix of the proposal distribution is \code{"SS'"}.
#'   \item{"updateLim":} Vector of two values: min and max length of updated state sequence
#'   \item{"updateProbs":} Probability for each element of \code{"updateLim[1]:updateLim[2]"} (if \code{"NULL"},
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
#' @param adapt Integer. (Experimental) If \code{"adapt"} > 0, use the the Robust Adaptive Metropolis (RAM) algorithm
#' by Vihola (2012) to update the proposal distribution for each parameter at each iteration (up to \code{"adapt"} iterations)
#' to a target acceptance rate of 23.4\%.
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
#'#load ctmm for data and comparison
#'library(ctmm)
#'#load pelican dataset
#'data('pelican')
#'#extract argos data with errors and store in a dataframe
#'Pelican <- data.frame(x = pelican$argos$x, y = pelican$argos$y, time = pelican$argos$t, ID = 1, error = pelican$argos$VAR.xy )
#'
#'#run 10000 iterations of a 2 state model (no error, but see commented line for example)
#'mcmc <- runMCMC( track = Pelican[,1:4], nbStates = 2, nbIter = 10000,
#'                 fixpar = list( tau_pos = c( NA, NA ), tau_vel = c( NA, NA ), sigma = c( NA, NA ) ),
#'                 inits = list( tau_pos = c( 13e6, 13e6 ),
#'                               tau_vel = c( 1e4, 1e4 ),
#'                               sigma = c( 2e11, 2e11 ),
#'                               Q = matrix( c( -0.05, 0.05, 0.05, -0.05 ), 2 ),
#'                               state = sample( 1:2, size = nrow( Pelican ), replace = TRUE ) ),
#'                 priors = list( mean = log( c( 9e6, 9e6, 1e4, 1e4, 2e11, 2e11 ) ),
#'                                sd = c( 2, 2, 2, 2, 2, 2 ), shape = 2, rate = 10, con = 0 ),
#'                 props = list( S = diag( c( 0.2, 0.2, 0.1, 0.1, 0.2, 0.2 ), 6 ), updateLim = c( 3, 150 ), updateProbs = rep( 1/148, 148 ) ),
#'                 tunes = list( thinStates = 10000 * 0.001 ),
#'                 #Hmat = cbind( Pelican$error, Pelican$error, matrix( rep( c( 0, 0 ), nrow( Pelican ) ), ncol = 2 ) ),
#'                 Hmat = matrix( rep( c( 0, 0, 0 ,0 ), nrow( Pelican ) ), ncol = 4 ) )
#'
#'#parameter estimates (tau_pos1, tau_pos2, tau_vel1, tau_vel2, sigma1, sigma2)
#'colMeans( mcmc$allparam[ -( 1:( nrow( mcmc$allparam ) / 2 ) ) , ] )
#'
#'#state sequence
#'round( colMeans( mcmc$allstates[ -( 1:nrow( mcmc$allstates ) / 2 ), ] ) )
#' }
#'
#' @importFrom stats dnorm runif rnorm rexp rgamma
#' @importFrom prettyunits vague_dt pretty_dt
#' @importFrom ramcmc adapt_S
#' @export
#'
#' @useDynLib MSctmm
runMCMC <- function( track, nbStates, nbIter, fixpar = NULL, inits, priors, props, tunes, Hmat,
                     updateState = TRUE, adapt = FALSE )
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
    param <- c( tau_pos, tau_vel, sigma )
    Q <- rep( list( inits$Q ), length( unique( track$ID ) ) )
    names(Q) <- unique( track$ID )
    state0 <- inits$state

    if( is.null( fixpar ) ) {
      fixpar <- list( tau_pos = rep( NA, nbStates ),
                      tau_vel = rep( NA, nbStates ),
                      sigma = rep( NA, nbStates ) )
    } else {
      fixpar <- fixpar
    }
    if(length( unlist( inits[1:3] ) )!=length( unlist( fixpar ) ) )
      stop( "'fixpar' has the wrong length" )

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
    S <- props$S
    updateLim <- props$updateLim
    updateProbs <- props$updateProbs
    if(is.null(S) | is.null(updateLim))
        stop("'props' should have components S and updateLim")

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

    ids <- as.character( unique( track$ID ) )
    obs <- list()
    switch <- list()
    data.list <- list()

    for( i in ids ) {
      nbObs <- nrow(track[ which( track$ID == i ), ])
      obs[[ i ]] <- matrix( c( track[ which( track$ID == i ), "x" ],
                               track[ which( track$ID == i ), "y" ],
                               track[ which( track$ID == i ), "ID" ],
                               track[ which( track$ID == i ), "time" ],
                               state0[ which( track$ID == i )  ] ),
                            ncol = 5 )
      colnames( obs[[ i ]] ) <- c( "x", "y", "ID", "time", "state" )
      indSwitch <- which( obs[[ i ]][ -1, "state" ] != obs[[ i ]][ -nbObs, "state" ] ) + 1
      switch[[ i ]] <- matrix( c( obs[[ i ]][ indSwitch, "time" ] - 0.001, rle( obs[[ i ]][ , "state" ] )$values[ -1 ] ), ncol = 2 )
      colnames( switch[[ i ]] ) <- c( "time", "state" )

      if( !all( is.na( switch[[ i ]] ) ) ) {
        data.list[[ i ]] <- rbind( obs[[ i ]] ,cbind( NA, NA, switch[[ i ]] ) )
      } else {
        data.list[[ i ]] <- obs[[ i ]]
      }
    }
    # flatten data
    data.df <- do.call( "rbind", data.list )
    data.df <- data.df[ order( data.df[,"ID"], data.df[,"time"] ), ]

    # initialise Hmat (rows of 0s for transitions)
    HmatAll <- matrix( 0, nrow( data.df ), 4 )
    HmatAll[ which( !is.na( data.df[ , "x" ] ) ), ] <- Hmat

    # initial likelihood
    oldllk <- kalman_rcpp( data = data.df, param = param, Hmat = HmatAll )[1]

    # initial log-prior
    oldlogprior <- sum( dnorm( log( param[ is.na(unlist(fixpar)) ] ), priorMean[ is.na(unlist(fixpar)) ], priorSD[ is.na(unlist(fixpar)) ], log = TRUE ) )

    ###############################
    ## Loop over MCMC iterations ##
    ###############################
    allparam <- matrix( NA, nrow = nbIter, ncol = 5 * nbStates )
    allrates <- array( NA, dim = c( nbIter, nbStates * ( nbStates - 1 ), length( ids ) ) )
    allstates <- matrix( NA, nrow = nbIter / thinStates, ncol = nrow( track ) ) # uses a lot of memory!
    accSwitch <- rep( 0, nbIter )
    accParam <- rep( 0, nbIter )
    allLen <- rep( NA, nbIter )

    t0 <- Sys.time()
    for( iter in 1:nbIter ) {
        if( iter %% 100 == 0 )
            cat( "\rIteration ", iter, "/", nbIter, "... ",
                vague_dt( difftime( Sys.time(), t0, units = "secs" ) / iter * ( nbIter - iter ), "short" ),
                " remaining (est)",
                " -- accSwitch = ", round( sum( accSwitch ) / iter * 100 ), "%",
                " -- accParam = ", round( sum( accParam ) / iter * 100 ), "%", sep = "" )

        ######################################
        ## 1. Update discrete state process ##
        ######################################
        if( updateState ) {

            # pick an individual at random for which we will update the state sequence
            id <- sample( ids, 1 )
            upState <- updateState( obs = obs[[ id ]], switch = switch[[ id ]], updateLim = updateLim,
                                   updateProbs = updateProbs, Q = Q[[ id ]])
            newData.list <- data.list
            newData.list[[ id ]] <- upState$newData
            newSwitch <- switch
            newSwitch[[ id ]] <- upState$newSwitch
            allLen[iter] <- upState$len

            # flatten data
            newData.df <- do.call( "rbind", newData.list )

            # update Hmat (rows of 0s for transitions)
            newHmatAll <- matrix( 0, nrow( newData.df ), 4 )
            newHmatAll[ which( !is.na( newData.df [ , "x" ] ) ), ] <- Hmat

            # Calculate acceptance ratio
            newllk <- kalman_rcpp( data = newData.df, param = param, Hmat = newHmatAll )[1]
            logHR <- newllk - oldllk

            if( log( runif(1) ) < logHR ) {
                # Accept new state sequence
                accSwitch[iter] <- 1
                switch <- newSwitch
                data.list <- newData.list
                obs <- lapply( data.list, function( data ) { data[!is.na( data[,"x"] ),] } )
                oldllk <- newllk
                HmatAll <- newHmatAll
            }

        }

        ###################################
        ## 2. Update movement parameters ##
        ###################################
        # On working scale
        u <- rnorm( length( param ) )
        thetas <- log( param ) + as.vector( S %*% u )
        newlogprior <- sum(dnorm(thetas[ is.na( unlist( fixpar ) ) ],priorMean[ is.na( unlist( fixpar ) ) ],priorSD[ is.na( unlist( fixpar ) ) ],log=TRUE))

        # On natural scale
        thetasprime <- unlist( fixpar )
        thetasprime[ is.na( unlist( fixpar ) ) ] <- exp( thetas[ is.na( unlist( fixpar ) ) ] )

        # Calculate acceptance ratio
        data.df <- do.call( "rbind", data.list )
        kalman <- kalman_rcpp( data = data.df, param = thetasprime, Hmat = HmatAll )
        newllk <- kalman[1]
        mu <- kalman[-1]
        logHR <- newllk + newlogprior - oldllk - oldlogprior

        if( log( runif( 1 ) ) < logHR ) {
            # Accept new parameter values
            accParam[iter] <- 1
            param <- thetasprime
            oldllk <- newllk
            oldlogprior <- newlogprior
        }

        if( adapt & iter >= 1000 & iter <= adapt ) {
            S[is.na(unlist(fixpar)), is.na(unlist(fixpar))] <- adapt_S(S[is.na(unlist(fixpar)), is.na(unlist(fixpar))], u[is.na(unlist(fixpar))], min( 1, exp(logHR) ), iter )
        }

        ###############################
        ## 3. Update switching rates ##
        ###############################
        Q <- lapply( ids, function( id ) { updateQ( nbStates = nbStates, data = data.list[[ id ]], switch = switch[[ id ]],
                     priorShape = priorShape, priorRate = priorRate,
                     priorCon = priorCon ) } )
        names(Q) <- ids

        #########################
        ## Save posterior draw ##
        #########################
        allparam[iter,] <- c( rbind( matrix( param, ncol = nbStates ), matrix( mu, ncol = nbStates ) ) )
        allrates[ iter, , ] <- matrix( unlist( lapply( Q, function( q ){ q[ !diag( nbStates ) ] } ) ), ncol = length( ids ), nrow = nbStates * ( nbStates - 1 ) )
        if( iter %% thinStates == 0 )
            allstates[iter / thinStates,] <- unlist( lapply( obs, function( ob ) { ob[ , "state" ] } ) )
    }
    cat( "\n" )
    cat( "Elapsed: ", pretty_dt( difftime( Sys.time(), t0, units = "secs" ) ), sep = "" )
    cat( "\n" )
    if( adapt ) {
      cat( 'updateLim: [', updateLim[1], ',', updateLim[2], ']\n', sep = '' )
      cat( 'S:\n')
      print( round( S, 4 ) )
      cat( '\n' )
    }

    return( list(allparam = allparam,
                allrates = allrates,
                allstates = allstates,
                accSwitch = accSwitch,
                allLen = allLen))
}
