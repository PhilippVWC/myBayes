#include <Rcpp.h>
using namespace Rcpp;

//' @title General symmetric map
//' @description a generalization of a one dimensional map, which incorporates (among others) the lorenz-, logistic and cubeMap
//' @param x input value for one iteration
//' @param r control parameter
//' @param alpha exponent
//' @return The result of one single iteration/
//' @author J.C. Lemm, P.v.W. Crommelin
//' @references S. Sprott, Chaos and Time-series analysis
//'

//' @export
// [[Rcpp::export]]
double gensymMap_cpp(double x, double r, double alpha){
  return(r*(1-pow(2*abs(x-0.5),alpha)));
}

//' @title Time series creation
//' @description Create time series produced by general symmetric map
//' @param N integer - Number of iterations
//' @param x0 double - starting value
//' @param r double - controll parameter
//' @param alpha double - exponent of general symmetric map
//' @param N_discr integer - controlls discretization of state space (zero corresponds to continuous case)
//' @param skipFirst Boolean - If set to FALSE, the resulting time series contains the initial value x0
//' @return vector of doubles - the resulting time series
//' @author J.C. Lemm, P. v.W. Crommelin
//' @export
// [[Rcpp::export]]
NumericVector gsm_iter_cpp(int N, double x0, double r, double alpha, int N_discr,bool skipFirst){
  Rcpp::NumericVector X(N);
  if(skipFirst){
    if(N_discr == 0){
      X[0] = gensymMap_cpp(x0,r,alpha);
    }else{
      X[0] = round(gensymMap_cpp(x0,r,alpha)*N_discr)/N_discr;
    }
  }
  if(N>1){
    if(N_discr == 0){
      for(int i = 1; i<N ; i++){
        X[i] = gensymMap_cpp(X[i-1],r,alpha);
      }
    }else{
      for(int i = 1; i<N ; i++){
        X[i] = round(gensymMap_cpp(X[i-1],r,alpha)*N_discr)/N_discr;
      }
    }
  }
  return(X);
}

/*** R
print("skipfirst true")
*/
