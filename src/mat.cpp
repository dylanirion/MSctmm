#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// dexp2
// numerically stable evaluation of 1 - exp(-x)^2
NumericVector dexp2( NumericVector x, Nullable<NumericVector> Exp = R_NilValue ) {
    NumericVector res( x.size() );
    for(int i=0; i<x.size(); i++) {
        if ( Exp.isNotNull() ) {
            NumericVector e(Exp);
            if( e(i) < 0.7071068 ) {
                res(i) = 1 - pow( e(i), 2 );
            } else {
                res(i) =  2 * e(i) * sinh( x(i) );
            }
        } else {
            NumericVector e = exp(-x);
            if( e(i) < 0.7071068 ) {
                res(i) = 1 - pow( e(i), 2 );
            } else {
                res(i) =  2 * e(i) * sinh( x(i) );
            }
        }
    }
    return res;
}

// dexp1
// numerically stable evaluation of 1 - exp(-x)^1
NumericVector dexp1( NumericVector x, Nullable<NumericVector> Exp = R_NilValue ) {
    NumericVector res(x.size());
    for(int i=0; i<x.size(); i++) {
        if (Exp.isNotNull()) {
            NumericVector e(Exp);
            if( e(i) < 0.5 ) {
                res(i) = 1 - e(i);
            } else {
                res(i) = 2 * sqrt( e(i) ) * sinh( x(i) / 2 );
            }
        } else {
            NumericVector e = exp(-x);
            if( e(i) < 0.5 ) {
                res(i) = 1 - e(i);
            } else {
                res(i) = 2 * sqrt( e(i) ) * sinh( x(i) / 2 );
            }
        }
    }
    return res;
}

// Sinc
// sinc function
NumericVector sinc( NumericVector x, Nullable<NumericVector> SIN = R_NilValue ) {
    NumericVector res(x.size());
    for(int i=0; i<x.size(); i++) {
        if (SIN.isNotNull()) {
            NumericVector s(SIN);
            if( x(i) == 0 ) {
                res(i) = 1;
            } else {
                res(i) = s(i) / x(i);
            }
        } else {
            NumericVector s = sin(x);
            if( x(i) == 0 ) {
                res(i) = 1;
            } else {
                res(i) = s(i) / x(i);
            }
        }
    }
    return res;
}

// Sinch
// hyperbolic sinc
NumericVector sinch( NumericVector x, Nullable<NumericVector> SINH = R_NilValue ) {
    NumericVector res(x.size());
    for(int i=0; i<x.size(); i++) {
        if (SINH.isNotNull()) {
            NumericVector s(SINH);
            if( x(i) == 0 ) {
                res(i) = 1;
            } else {
                res(i) = s(i) / x(i);
            }
        } else {
            NumericVector s = sinh(x);
            if( x(i) == 0 ) {
                res(i) = 1;
            } else {
                res(i) = s(i) / x(i);
            }
        }
    }
    return res;
}
// Make T matrix
//
// This code is adapted from the langevin function in the R package ctmm (Calabrese et al., 2016).
mat makeT(double tau_pos, double tau_vel, double dt) {

    vec tau(2);
    tau(0) = tau_pos;
    tau(1) = tau_vel;

    double omega2;
    double f;
    double nu;
    vec dtau;
    vec Exp;
    vec DExp;

    //IID limit
    mat T(4,4, fill::zeros);

    if(std::isinf(tau(0)) && tau(1) == 0) { // BM
        if(!std::isinf(dt)) {
          T(0,0) = 1;
        }
    } else if(!std::isinf(tau(0)) && tau(1) == 0) { // OU
      if(!std::isinf(dt)) {
        dtau = dt / tau(0);
        Exp = exp( -dtau );
        T(0,0) = as<double>(wrap(Exp));
      }
    } else if(std::isinf(tau(0)) && tau(1) != 0) { // IOU
      omega2 = 1/tau(1);
      f = as_scalar( mean( 1/tau ) );
      nu = as_scalar( diff( 1 / tau ) ) / 2;
      dtau = dt / tau(1);
      Exp = exp( -dtau );
      DExp = as<vec>( dexp1( as<NumericVector>( wrap( dtau ) ), as<NumericVector>( wrap ( Exp ) ) ) );

      if(!std::isinf(dt)) {
        T(0,0) = 1;
        T(0,1) = tau(1) * as_scalar( DExp );
        T(1,1) = as_scalar( Exp );
      }
    } else if( tau(0) > tau(1) ){ // overdamped
      f = as_scalar( mean( 1/tau ) );
      nu = as_scalar( diff( 1 / tau ) ) / 2;
      omega2 = ( 1 / tau(0) ) * ( 1 / tau(1) );
    } else if( tau(0) == tau(1) ) { // critically damped
      f = 1/tau(0);
      nu = 0;
      omega2 = ( 1 / tau(0) ) * ( 1 / tau(1) );
    }

    double fdt = f * dt;

    if(!std::isinf(tau(0)) && tau(0) != 0 && tau(1) != 0) { //IOU, OUF/OUO
      if(!std::isinf(dt)) {

        double nudt = nu * dt;
        //bool EXP = tau(0) > tau(1) && nudt > 0.8813736;
        bool EXP = nudt > 0.8813736;

        double c0;
        double c1;
        double c2;

        if( EXP ) { // Exponential functions
          dtau = dt / tau;
          double dift = as_scalar( diff( tau ) );
          vec Exp0 = exp( -dtau );
          Exp = Exp0 / dift;
          c0 = as_scalar( diff( Exp % tau ) );
          c1 = as_scalar( -diff( Exp ) );
          c2 = as_scalar( diff( Exp / tau ) );
        } else { // Trigonometric and hyperbolic-trigonometric functions
          Exp = exp( -fdt );
          double Sin0;
          double Sinc0;
          double Cos0;

          //if( tau(0) > tau(1) ) { // Hyperbolic-trigonometric
          Sin0 = sinh( nudt );
          Sinc0 = as<double>( wrap( sinch( as<NumericVector>( wrap( nudt ) ), as<NumericVector>( wrap ( Sin0 ) ) ) ) );
          Cos0 = cosh( nudt );
          /*
          } else { // Trigonometric
           Sin0 = sin(nudt);
           Sinc0 = as<double>( wrap( sinc( as<NumericVector>( wrap( nudt ) ), as<NumericVector>( wrap( Sin0 ) ) ) ) );
           Cos0 = cos(nudt);
          }
           */
          double SincE = Sinc0 * as_scalar( Exp );
          double CosE = Cos0 * as_scalar( Exp );

          c0 = CosE + fdt * SincE;
          c1 = -( omega2 * dt ) * SincE;
          c2 = -omega2 * ( CosE - fdt * SincE );
        }

        T(0,0) = c0;
        T(1,0) = c1;
        T(0,1) = -c1 / omega2;
        T(1,1) = -c2 / omega2;
      }
    }

    T.submat(2,2,3,3) = T.submat(0,0,1,1);

    return T;
}

// Make Q matrix
//
// This code is adapted from the langevin function in the R package ctmm (Calabrese et al., 2016).
mat makeQ(double tau_pos, double tau_vel, arma::mat sigma, double dt) {

    vec tau(2);
    tau(0) = tau_pos;
    tau(1) = tau_vel;

    double omega2;
    double f;
    double nu;
    double tfomega2;
    vec dtau;
    vec Exp;
    vec DExp;
    vec DExp2;

    //IID limit
    mat Q(4,4, fill::zeros);
    Q(0,0) = 1;

    if(std::isinf(tau(0)) && tau(1) == 0) { // BM
      // absorbing 1/tau into sigma # VAR -> Diffusion
      Q(0,0) = 2 * dt;
    } else if(!std::isinf(tau(0)) && tau(1) == 0) { // OU
      if(!std::isinf(dt)) {
        dtau = dt / tau(0);
        Exp = exp( -dtau );
        Q(0,0) = as<double>(wrap(dexp2(as<NumericVector>(wrap(dtau)), as<NumericVector>(wrap(Exp)))));
        omega2 = 1/tau(0);
        Q(1,1) = omega2;
      }
    } else if(std::isinf(tau(0))) { // IOU
      omega2 = 1/tau(1);
      Q(1,1) = omega2;
      f = as_scalar(1 / tau(1) / 2);
      nu = as_scalar(1 / tau(1) / 2);
      tfomega2 = 2 * f / omega2;
    } else if( tau(0) > tau(1) ){ // overdamped
      f = as_scalar( mean( 1/tau ) );
      nu = as_scalar( diff( 1 / tau ) ) / 2;
      omega2 = ( 1 / tau(0) ) * ( 1 / tau(1) );
      tfomega2 = 2 * f / omega2;
    } else if( tau(0) == tau(1) ) { // critically damped
      f = 1/tau(0);
      nu = 0;
      omega2 = ( 1 / tau(0) ) * ( 1 / tau(1) );
      tfomega2 = 2 * f / omega2;
    }

    double fdt = f * dt;

    if(std::isinf(tau(0)) && tau(1) != 0) { //IOU
      dtau = dt / tau(1);
      Exp = exp( -dtau );
      DExp = as<vec>( wrap( dexp1( as<NumericVector>( wrap( dtau ) ), as<NumericVector>( wrap( Exp ) ) ) ) );
      DExp2 = pow( DExp, 2 ); // (1-exp(-dt/tau[2]))^2
      Q(0,0) = as_scalar( clamp( 2 * dt - tau(1) * ( 2 * DExp + DExp2 ), 0.0, datum::inf ) );
      Q(1,1) = as<double>( wrap( dexp2( as<NumericVector>( wrap( dtau ) ), as<NumericVector>( wrap ( Exp ) ) ) ) ) / tau(1);

      if( !std::isinf( dt ) ) {
        Q(1,0) = as_scalar( clamp( DExp2, 0.0, sqrt( Q(0,0) * Q(1,1) ) ) );
        Q(0,1) = as_scalar( clamp( DExp2, 0.0, sqrt( Q(0,0) * Q(1,1) ) ) );
      }
    } else if(tau(0) != 0 && tau(1) != 0) {
      if( !std::isinf( dt ) ) { //IOU,OUF/OUO,IID

        double nudt = nu * dt;
        //bool EXP = tau(0) > tau(1) && nudt > 0.8813736;
        bool EXP = nudt > 0.8813736;

        double c1;

        // Initially cancelling terms
        if( EXP ) { // Exponential functions
          dtau = dt / tau;
          double dift = as_scalar( diff( tau ) );
          vec Exp0 = exp( -dtau );
          Exp = Exp0 / dift;
          vec T2 = pow( tau, 2 );
          c1 = -as_scalar( diff( Exp ) );
          double dift2 = pow( dift, 2 );
          double S1 = as<double>( wrap( dexp2( as<NumericVector>( wrap( dtau(0) ) ), as<NumericVector>( wrap( Exp0(0) ) ) ) ) );
          double S2 = as<double>( wrap( dexp2( as<NumericVector>( wrap( dtau(1) ) ), as<NumericVector>( wrap( Exp0(1) ) ) ) ) );
          double S12 = 2 * tau(0) * tau(1) * as<double>( wrap( dexp1( as<NumericVector>( wrap( fdt ) ), as<NumericVector>( wrap( Exp0(0) * Exp0(1) ) ) ) ) );

          Q(0,0) = ( T2(0) * S1 - S12 + T2(1) * S2) / dift2;
          Q(1,1) = ( T2(1) * S1 - S12 + T2(0) * S2) / dift2 * omega2;
        } else { // Trigonometric and hyperbolic-trigonometric functions
          //Need these again
          Exp = exp( -fdt );
          double Sin0;
          double Sinc0;
          double Cos0;

          //if( tau(0) > tau(1) ) { // Hyperbolic-trigonometric
          Sin0 = sinh( nudt );
          Sinc0 = as<double>( wrap( sinch( as<NumericVector>( wrap ( nudt ) ), as<NumericVector>( wrap( Sin0 ) ) ) ) );
          Cos0 = cosh( nudt );
          /*
          } else { // Trigonometric
           Sin0 = sin( nudt );
           Sinc0 = as<double>( wrap( sinc( as<NumericVector>( wrap ( nudt ) ), as<NumericVector>( wrap( Sin0 ) ) ) ) );
           Cos0 = cos( nudt );
          }
           */
          // Also need these again
          double SincE = Sinc0 * as_scalar( Exp );
          c1 = -( omega2 * dt ) * SincE;

          double CROSS = fdt * Sinc0 * as_scalar(  Exp );
          double OUTER = pow( Cos0, 2 ) * as<double>( wrap( dexp2( as<NumericVector>( wrap( fdt ) ), as<NumericVector>( wrap( Exp ) ) ) ) ) - pow( CROSS, 2 );
          CROSS = 2 * Cos0 * as_scalar( Exp ) * CROSS;
          double Sin2 = pow( Sin0, 2 );

          //if( tau(0) > tau(1) )
          //{
            Q(0,0) = OUTER - Sin2 - CROSS;
            Q(1,1) = ( OUTER - Sin2 + CROSS ) * omega2;
          //}
          //else
          //{
          //  Q(0,0) = OUTER + Sin2 - CROSS;
          //  Q(1,1) = ( OUTER + Sin2 + CROSS ) * omega2;
          //}
        }

        //Initially vanishing terms
        double c12 = pow( c1, 2 );
        double TT = tfomega2;
        Q(0,0) = ( Q(0,0) - c12 / omega2 );
        Q(0,1) = TT * c12;
        Q(1,0) = TT * c12;
        Q(1,1) = ( Q(1,1) - c12 );

      } // end OUF/OUO
    }
    Q.submat(0,0,1,1) = Q.submat(0,0,1,1) * sigma.t();

    //IOU prior fix
    Q.replace(datum::nan, 0 );  // replace each NaN with 0
    Q.submat(2,2,3,3) = Q.submat(0,0,1,1);

    return Q;
}
