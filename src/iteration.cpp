#include <Rcpp.h>
using namespace Rcpp;

//' @title General symmetric map
//' @name gensymMap_cpp
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

//' @export
// [[Rcpp::export]]
NumericVector gsm_iter_cpp(int N, double x0, double r, double alpha, bool N_discr,bool skipFirst){
  printf("printing");
  NumericVector X(N);
  if(skipFirst){
 /*** R
 print("skipfirst true")
 */
    printf("skipfirst true");
    if(N_discr == NULL);
  }
}
