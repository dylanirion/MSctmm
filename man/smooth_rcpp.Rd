% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{smooth_rcpp}
\alias{smooth_rcpp}
\title{Kalman smoother}
\usage{
smooth_rcpp(data, nbStates, param, fixmu, Hmat)
}
\arguments{
\item{data}{Matrix of data, including columns \code{x}, \code{y}, \code{time}, \code{ID} and \code{state} (in that order).}

\item{param}{Vector of movement parameters (\code{tau_vel}, \code{tau_pos}, and \code{sigma})}

\item{fixmu}{Vector of mean locations for the OUF process (\code{x}, \code{y})}

\item{Hmat}{Matrix of observation error variance (four columns, and one row
for each row of data)}
}
\value{
a named List containing the predicted locations and velocities, and variance of these estimates.
}
\description{
This code is adapted from the package ctmm (Calabrese et al., 2016) crawl (Johnson et al., 2008),
and MScrawl (Michelot and Blackwell, 2019).
}
\references{
Calabrese, J.M., Fleming, C.H. and Gurarie, E. (2016).
ctmm: an r package for analyzing animal relocation data as a continuous‐time stochastic process.
Methods Ecol Evol, 7: 1124-1132. doi:10.1111/2041-210X.12559

Fleming, C.H., Sheldon, D., Gurarie, E., Fagan, W.F., LaPoint, S., Calabrese, J.M. (2017).
Kálmán filters for continuous-time movement models.
Ecol Inform, 40: 8-21. doi:10.1016/j.ecoinf.2017.04.008

Johnson, D.S., London, J.M., Lea, M.A., and Durban, J.W. (2008).
Continuous-time correlated random walk model for animal telemetry data.
Ecology, 89: 1208-1215. doi:10.1890/07-1032.1

Michelot, T., Blackwell, P.G. (2019).
State‐switching continuous‐time correlated random walks.
Methods Ecol Evol, 10: 637-649. doi:10.1111/2041-210X.13154
}
