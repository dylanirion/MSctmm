
#' Run MCMC iterations
#'
#' @param track Dataframe of data, with columns \code{"x"}, \code{"y"}, \code{"time"}, and \code{"ID"}
#' @param nbStates Number of states
#' @param nbIter Number of iterations
#' @param fixPar List of fixed parameter values (\code{"tau_pos"}, \code{"tau_vel"}, \code{"sigma"}), with \code{"NA"} for parameters to estimate
#' @param fixMu List of fixed OUF range centre coordinate pairs, with \code{"NA"} for IOU states or pairs to estimate
#' @param inits List of initial parameters (\code{"tau_pos"}, \code{"tau_vel"}, \code{"sigma"}, \code{"Q"}, \code{"state"})
#' @param priors List of parameters of prior distributions, with components:
#' \itemize{
#'   \item{"mean":} Vector of means for normal priors on movement parameters, of length \code{"3*nbStates"}
#'   \item{"sd":} Vector of standard deviations for normal priors on movement parameters, of length \code{"3*nbStates"}
#'   \item{"shape":} Vector of shapes of gamma priors for the transition rates
#'   \item{"rate":} Vector of rates of gamma priors for the transition rates
#'   \item{"con":} Vector of concentrations of Dirichlet priors for transition probabilities
#' }
#' @param knownStates Vector of known states with \code{"NA"} or those to estimate.
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
#' @param updateState Logical. If \code{"FALSE"}, the state process is not updated
#' (for exploratory analysis only, useful for testing single state models)
#' @param adapt Integer. If \code{"adapt"} > 0, use the the Robust Adaptive Metropolis (RAM) algorithm
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
#' }
#'
#' @importFrom stats dnorm runif rnorm rexp rgamma
#' @importFrom prettyunits vague_dt pretty_dt
#' @importFrom ramcmc adapt_S
#' @export
#'
#' @useDynLib MSctmm
runMCMC <- function( track, nbStates, nbIter, fixPar = NULL, fixMu = NULL, inits, priors,
                     knownStates, props, tunes, Hmat,
                     updateState = TRUE, adapt = FALSE )
{

    ######################
    ## Unpack arguments ##
    ######################
    # initial parameters
    tau_pos <- inits$tau_pos
    tau_vel <- inits$tau_vel
    #NB considering only isotropic case for now, this will eventually need to be a matrix
    #kalman_rcpp and MakeQ/Sigma will need to accept a matrix
    sigma <- inits$sigma
    param <- c( tau_pos, tau_vel, sigma )
    state0 <- inits$state
    #overwrite state inits with knownStates
    state0[which( !is.na( knownStates ) )] <- knownStates[which( !is.na( knownStates ) )]


    if( is.null( fixPar ) ) {
      fixpar <- list( tau_pos = rep( NA, nbStates ),
                      tau_vel = rep( NA, nbStates ),
                      sigma = rep( NA, nbStates ) )
    } else {
      fixpar <- fixPar
    }
    if( length( unlist( inits[1:3] ) ) != length( unlist( fixpar ) ) )
      stop( "'fixpar' has the wrong length" )

    if( is.null( fixMu ) ) {
      fixmu <- rep( list( c( NA, NA) ), nbStates )
    } else {
      fixmu <- fixMu
    }
    if( 2 * nbStates != length( unlist( fixmu ) ) )
      stop( "'fixMu' has the wrong length" )

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
    thinStates <- ifelse( is.null( tunes$thinStates ), 1, tunes$thinStates )

    ####################
    ## Initialisation ##
    ####################
    # Prepare data structures
    if( !( "ID" %in% colnames( track ) ) ) {
      warning( "'track' should have column 'ID', assuming data is one individual"  )
      track$ID <- 1
    }

    track$ID <- as.numeric( as.factor( track$ID ) )

    ids <- unique( track$ID )
    obs <- list()
    known <- list()
    switch <- list()
    data.list <- list()

    if( length(inits$Q) == length(ids) & "list" %in% class(inits$Q) ) {
      Q <- inits$Q
    } else {
      Q <- rep( list( inits$Q ), length( unique( track$ID ) ) )
    }
    names(Q) <- ids

    if(is.null(tau_pos) | is.null(tau_vel) | is.null(sigma) | is.null(Q) | is.null(state0) )
      stop("'inits' should have components tau_pos, tau_vel, sigma, Q, and state")

    for( id in ids ) {
      nbObs <- nrow(track[ which( track$ID == id ), ])
      obs[[ id ]] <- matrix( c( track[ which( track$ID == id ), "x" ],
                               track[ which( track$ID == id ), "y" ],
                               track[ which( track$ID == id ), "time" ],
                               as.numeric( track[ which( track$ID == id ), "ID" ] ),
                               state0[ which( track$ID == id )  ] ),
                            ncol = 5 )
      known[[ id ]] <- knownStates[ which( track$ID == id )  ]
      colnames( obs[[ id ]] ) <- c( "x", "y", "time", "ID",  "state" )
      indSwitch <- which( obs[[ id ]][ -1, "state" ] != obs[[ id ]][ -nbObs, "state" ] ) + 1
      switch[[ id ]] <- matrix( c( obs[[ id ]][ indSwitch, "time" ] - 0.001, rle( obs[[ id ]][ , "state" ] )$values[ -1 ] ), ncol = 2 )
      colnames( switch[[ id ]] ) <- c( "time", "state" )

      if( !all( is.na( switch[[ id ]] ) ) ) {
        data.list[[ id ]] <- rbind( obs[[ id ]][,c( "x", "y", "time", "ID",  "state" )] ,cbind( NA, NA, switch[[ id ]][,"time"], id, switch[[ id ]][,"state"] ) )
        data.list[[ id ]] <- data.list[[ id ]][ order( data.list[[ id ]][,"time"] ),]
      } else {
        data.list[[ id ]] <- obs[[ id ]][,c( "x", "y", "time", "ID",  "state" )]
      }
    }
    # flatten data
    data.mat <- do.call( "rbind", data.list )
    data.mat <- data.mat[ , c( "x", "y", "time", "ID", "state" ) ]

    # initialise Hmat (rows of 0s for transitions)
    HmatAll <- matrix( 0, nrow( data.mat ), 4 )
    HmatAll[ which( !is.na( data.mat[ , "x" ] ) ), ] <- Hmat

    # initial likelihood
    oldllk <- kalman_rcpp( data = data.mat, param = param, fixmu = unlist( fixmu ), Hmat = HmatAll )[1]

    # initial log-prior
    oldlogprior <- sum( dnorm( log( param[ is.na( unlist( fixpar ) ) ] ), priorMean[ is.na( unlist( fixpar ) ) ], priorSD[ is.na( unlist( fixpar ) ) ], log = TRUE ) )

    ###############################
    ## Loop over MCMC iterations ##
    ###############################
    allparam <- matrix( NA, nrow = nbIter, ncol = 5 * nbStates )
    colnames(allparam) <- c( paste("tau_pos[", 1:nbStates, "]", sep = ""),
                             paste("tau_vel[", 1:nbStates, "]", sep = ""),
                             paste("sigma[", 1:nbStates, "]", sep = ""),
                             paste( c( "mu_x[", "mu_y["), rep(1:nbStates, each = 2 ), c( "]", "]" ), sep = "") )
    accParam <- rep( 0, nbIter )
    allrates <- array( NA, dim = c( nbIter, nbStates * ( nbStates - 1 ), length( ids ) ) )
    allstates <- matrix( NA, nrow = nbIter / thinStates, ncol = nrow( track ) ) # uses a lot of memory!
    accSwitch <- rep( 0, nbIter )
    allLen <- matrix( NA, nrow = nbIter, ncol = length( ids ) )
    allnLLk <- rep( NA, nbIter )

    t0 <- Sys.time()
    for( iter in 1:nbIter ) {
        if( iter %% 100 == 0 ) {
            cat( "\rIteration ", iter, "/", nbIter, "... ",
                vague_dt( difftime( Sys.time(), t0, units = "secs" ) / iter * ( nbIter - iter ), "short" ),
                " remaining (est)",
                " -- accSwitch = ", round( sum( accSwitch ) / iter * 100 ), "%",
                " -- accParam = ", round( sum( accParam ) / iter * 100 ), "%",
                "          ", sep = "" )
        }
        ######################################
        ## 1. Update discrete state process ##
        ######################################
        if( updateState ) {
          for( id in ids ) {
            upState <- updateState( obs = obs[[ id ]], knownStates = known[[ id ]], switch = switch[[ id ]], updateLim = updateLim,
                                    updateProbs = updateProbs, Q = Q[[ id ]])
            newData.list <- data.list
            newData.list[[ id ]] <- upState$newData[,c( "x", "y", "time", "ID",  "state" )]
            newData.list[[ id ]] <- newData.list[[ id ]][ order( newData.list[[id]][,"time"] ),]
            newSwitch <- switch
            newSwitch[[ id ]] <- upState$newSwitch
            allLen[iter, which(ids == id) ] <- upState$len
          }

          # flatten data
          newData.mat <- do.call( "rbind", newData.list )
          newData.mat <- newData.mat[ , c( "x", "y", "time", "ID", "state" ) ]

          # update Hmat (rows of 0s for transitions)
          newHmatAll <- matrix( 0, nrow( newData.mat ), 4 )
          newHmatAll[ which( !is.na( newData.mat[ , "x" ] ) ), ] <- Hmat

          # Calculate acceptance ratio
          newllk <- kalman_rcpp( data = newData.mat, param = param, fixmu = unlist( fixmu ), Hmat = newHmatAll )[1]
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
        pass = F
        while( !pass ) {  # ensure tau_p >= tau_v
          # On working scale [-Inf,Inf]
          u <- rnorm( length( param ) )
          thetas <- log( param ) + as.vector( S %*% u )
          # On natural scale [0, Inf ]
          thetasprime <- unlist( fixpar )
          thetasprime[ is.na( unlist( fixpar ) ) ] <- exp( thetas[ is.na( unlist( fixpar ) ) ] )

          # is there potential to get stuck here?
          if( all( thetasprime[1:nbStates] >= thetasprime[( nbStates + 1 ):( 2 * nbStates )] ) ) {
            pass <- T
          }
        }
        newlogprior <- sum( dnorm( thetas[ is.na( unlist( fixpar ) ) ], priorMean[ is.na( unlist( fixpar ) ) ], priorSD[ is.na( unlist( fixpar ) ) ], log = TRUE ) )

        # Calculate acceptance ratio
        data.mat <- do.call( "rbind", data.list )
        data.mat <- data.mat[ , c( "x", "y", "time", "ID", "state" ) ]
        kalman <- kalman_rcpp( data = data.mat, param = thetasprime, fixmu = unlist( fixmu ), Hmat = HmatAll )
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
          #S[is.na(unlist(fixpar)), is.na(unlist(fixpar))] <- adapt_S(S[is.na(unlist(fixpar)), is.na(unlist(fixpar))], u[is.na(unlist(fixpar))], min( 1, exp(logHR) ), iter )
          # calculate S by state instead
          for( i in 1:nbStates ) {
            index <- seq( i, length( param ), by = nbStates )
            S[index,index][is.na(unlist(fixpar)[index]), is.na(unlist(fixpar)[index])] <- adapt_S(S[index,index][is.na(unlist(fixpar)[index]), is.na(unlist(fixpar)[index])], u[index][is.na(unlist(fixpar)[index])], min( 1, exp(logHR) ), iter )
          }
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
        allparam[iter,] <- cbind( matrix( param, ncol = 3 * nbStates ), matrix( mu, ncol =  2 * nbStates ) )
        allrates[iter, , ] <- matrix( unlist( lapply( Q, function( q ){ q[ !diag( nbStates ) ] } ) ), ncol = length( ids ), nrow = nbStates * ( nbStates - 1 ) )
        if( iter %% thinStates == 0 ){
          allstates[iter / thinStates,] <- unlist( lapply( obs, function( ob ) { ob[ , "state" ] } ) )
        }
        allnLLk[iter] <- oldllk
    }
    cat( "\n" )
    cat( "Elapsed: ", pretty_dt( difftime( Sys.time(), t0, units = "secs" ) ), sep = "" )
    cat( "\n" )

    return( list( allparam = allparam,
                  allrates = allrates,
                  allstates = allstates,
                  accSwitch = accSwitch,
                  allLen = allLen ) )
}
