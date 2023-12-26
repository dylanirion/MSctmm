#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Adjoint of matrix
mat Adj(mat M){
  return conj(M).t();
}

// Hermitian part of matrix
mat He(mat M){
  return (M + Adj(M))/2;
}

// cpp implementation of Positive definite solver
// https://github.com/ctmm-initiative/ctmm/blob/73e7c00179eccb3ba888eb96d31cfe85baf111c4/R/matrix.R#L256
mat PDsolve(mat M, bool sym=true, bool force=false, bool pseudo=false, double tol=std::numeric_limits<double>::epsilon()){

  uvec DIM = {M.n_rows, M.n_cols};

  // check for Inf & invert those to 0 (and vice versa)
  bool ANYINF = M.diag().has_inf();
  bool ANYZERO = any(M.diag() <= 0) && sym;
  if(ANYINF || ANYZERO){
    // 1/Inf == 0 # correlations not accounted for
    if(ANYINF){
      M.each_row(find_nonfinite(M.diag())) = zeros(1, DIM(1));
      M.each_col(find_nonfinite(M.diag())) = zeros(DIM(0), 1);
    }

    // 1/0 == Inf
    if(ANYZERO){
      uvec i = find(M.diag() <= 0);
      M.each_row(i) = zeros(1, DIM(1));
      M.each_col(i) = zeros(DIM(0), 1);
      for(uword j = 0; j < i.n_elem; j++){
        M(i(j), i(j)) = datum::inf;
      }
    }

    // regular inverse of remaining dimensions
    std::vector<int> vREM;
    uvec finite = find_finite(M.diag());
    uvec zero = find(M.diag() <= 0);
    std::sort(finite.begin(), finite.end());
    std::sort(zero.begin(), zero.end());
    std::set_difference(finite.begin(), finite.end(), zero.begin(), zero.end(), std::back_inserter(vREM));
    uvec REM = arma::conv_to<uvec>::from(vREM);
    if(REM.size() != 0 ){
      M.submat(REM,REM) = PDsolve(M.submat(REM,REM),force=force, pseudo=pseudo, sym=sym);
    }

    return M;
  }
  if(!force && !pseudo){
    if(DIM(0) == 1){
      M = 1/M;

      return M;
    }

    if(DIM(0) == 2){
      double DET = M(0,0)*M(1,1)-M(0,1)*M(1,0);
      if(DET <= 0){ return diagmat(vec(2).fill(datum::inf)); } // force positive definite / diagonal
      double SWP = M(0,0);
      M(0,0) = M(1,1);
      M(1,1) = SWP;
      M(0,1) = -M(0,1);
      M(1,0) = -M(1,0);
      M = M/DET;

      return M;
    }
  }

  // symmetrize
  if(sym){ M = He(M); }

  // rescale
  vec W = abs(M.diag());
  W = sqrt(W);
  uvec zero = find(W <= tol);
  uvec notzero = find(W > tol);
  if(zero.size() != 0){ // do not divide by zero or near zero
    if(zero.size() != 0){
      for(uword j = 0; j < zero.n_elem; j++){
        W(zero(j)) = min(W.elem(notzero));  // do some rescaling... assuming axes are similar
      }
    }
    else{
      for(uword j = 0; j < zero.n_elem; j++){
        W(zero(j)) = 1;  // do no rescaling
      }
    }
  }

  mat W2 = W * W.t();

  // now a correlation matrix that is easier to invert
  M = M/W2;

  // try ordinary inverse
  M = solve(M, eye(DIM(0), DIM(1)));

  // back to covariance matrix
  M = M/W2;

  // symmetrize
  if(sym) { M = He(M); }

  return M;
}
