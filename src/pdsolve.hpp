#ifndef _PDSOLVE_
#define _PDSOLVE_

arma::mat PDsolve(arma::mat M, bool sym=true, bool force=false, bool pseudo=false, double tol=std::numeric_limits<double>::epsilon());

#endif
