#' Simulate from 2D IOU/OUF process
#'
#' @param obsTimes Vector of observation times
#' @param tau_pos Parameter \eqn{tau[pos]} of the movement process
#' @param tau_vel Parameter \eqn{tau[vel]} of the movement process
#' @param sigma Parameter \eqn{sigma} of the movement process
#' @param Q Infinitesimal generator matrix
#'
#' @return Simulated track
#'
#' @export
sim <- function( obsTimes, tau_pos, tau_vel, sigma, Q ) {
  nbStates <- length( tau_pos )
  nbObs <- length( obsTimes )
  obsTimes <- obsTimes - obsTimes[1]

  ############################
  ## Simulate state process ##
  ############################
  S <- NULL
  switchTimes <- NULL
  st <- 1 # start in state 1
  if ( nbStates > 1 ) {
    t <- rexp( 1, -Q[st,st] )
    while( t < obsTimes[nbObs] ) {
      if( nbStates > 2 ) {
        probs <- Q[st,-st] / sum( Q[st,-st] )
        st <- sample( ( 1:nbStates )[-st], size = 1, prob = probs )
      } else {
        st <- (1:2)[-st]
      }

      S <- c( S, st )
      switchTimes <- c( switchTimes, t )
      t <- t + rexp( 1, -Q[st,st] )
    }
  }

  times <- c( obsTimes, switchTimes )
  nbData <- length(times)

  S <- c( rep( NA, nbObs ), S )
  S <- S[order(times)]
  S[1] <- 1

  for( t in 1:nbData )
    if( is.na( S[t] ) )
      S[t] <- S[t-1]

  times <- sort(times)
  isObs <- ( times %in% obsTimes )
  dt <- c( Inf, diff( times ) )
  nbData <- length(times)

  ##############################################
  ## Simulate velocity and position processes ##
  ##############################################
  data <- matrix( 0, nrow = nbData, ncol = 4 )
  H <- c( 0, 0 , 0, 0 )
  for( t in 1:nbData )  {
    s <- S[t]
    Sigma <- makeSigma( tau_pos[s], tau_vel[s], sigma[s], dt[t] )
    if( t == 1 & is.infinite( tau_pos[s] ) ) {
      Sigma[] <- 0
    }
    Sigma <- PDfunc( Sigma, func = function(x){ sqrt( abs( x ) ) }, pseudo = TRUE )
    H <- makeMu( tau_pos[s], tau_vel[s], dt[t] ) %*% H + Sigma %*% rnorm( 4 )
    data[t,] <- H
  }

  return( data.frame( x = data[,1], y = data[,3], time = times, state = S, vx = data[,2], vy = data[,4] )[isObs,] )
}

#### borrowed CTMM functions
# map function for real-valued PSD matrices
#' @author Christen H Fleming, \email{flemingc@@si.edu}
PDfunc <-function(M,func=function(m){1/m},force=FALSE,pseudo=FALSE,tol=.Machine$double.eps)
{
  DIM <- dim(M)[1]
  if(is.null(DIM)) { DIM <- 1 }
  # tol <- max(tol,.Machine$double.eps)

  if(DIM==1)
  { M <- c(M) }
  else if(DIM==2)
  {
    TR <- (M[1,1] + M[2,2])/2 # half trace
    BIGNUM <- TR^2 > .Machine$double.xmax * .Machine$double.eps

    if(BIGNUM)
    {
      DET <- (M[1,1]/TR)*(M[2,2]/TR) - (M[1,2]/TR)*(M[2,1]/TR)
      DET <- 1 - DET # (tr^2 - det)/tr^2
    }
    else
    {
      DET <- M[1,1]*M[2,2] - M[1,2]*M[2,1]
      DET <- TR^2 - DET # tr^2 - det
    }

    if(DET<=0) # det is too close to tr^2
    {
      M <- diag(M)
      V <- array(0,c(2,2,2))
      V[1,1,1] <- V[2,2,2] <- 1
    }
    else
    {
      DET <- sqrt(DET) # now root term
      if(BIGNUM) { DET <- DET * TR }

      # hermitian formula
      V <- diag(1/2,2) %o% c(1,1) + ((M-TR*diag(2))/(DET*2)) %o% c(1,-1)

      M <- TR + c(1,-1)*DET
    }
  }
  else if(DIM>2) # arbitrary DIM
  {
    M <- eigen(M)
    V <- M$vectors
    M <- Re(M$values)

    V <- vapply(1:DIM,function(i){Re(V[,i] %o% Conj(V[,i]))},diag(1,DIM))
  }

  if(any(M<0) && !force && !pseudo) { stop("Matrix not positive definite.") }

  # negative eigenvalues indicate size of numerical error
  # MIN <- last(M)
  # if(MIN<0) { tol <- max(tol,2*abs(MIN)) }

  FORCE <- (M < tol) -> PSEUDO
  # PSEUDO <- (abs(M) < tol) # why abs(M)?

  if(any(FORCE) && force) { M[FORCE] <- tol }
  M <- func(M)
  if(any(PSEUDO) && pseudo) { M[PSEUDO] <- 0 }

  if(DIM==1)
  { M <- cbind(M) }
  else
  {
    # add up from smallest contribution to largest contribution
    INDEX <- sort(abs(M),method="quick",index.return=TRUE)$ix
    M <- lapply(INDEX,function(i){nant(M[i]*V[,,i],0)}) # 0/0 -> 0
    M <- Reduce('+',M)
  }

  return(M)
}

# 0/0 -> NaN -> to
# fixes a priori known limits
nant <- function(x,to)
{
  NAN <- is.nan(x)
  if(any(NAN)) { x[NAN] <- to }
  return(x)
}
