#include <Rcpp.h>
using namespace Rcpp;


//' @title Running average
//' @description This function calculates the running average of a given data vector
//' @param v Vector of type double - data vector
//' @return Vector of type double - The running mean/average
//' @author Philipp van Wickevoort Crommelin
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ravrg(Rcpp::NumericVector v){
  int N = v.size();
  Rcpp::NumericVector cumv(N);
  Rcpp::NumericVector cumvNormed(N);
  cumv[0] = v[0];
  cumvNormed[0] = cumv[0];
  for( int i=1 ; i<N ; i++ ){
    cumv[i] = v[i] + cumv[i-1];
    cumvNormed[i] = cumv[i]/(i+1);
  }
  return(cumvNormed);
}


//' @title Get index position on discrete grid
//' @description This function provides the corresponding index position on the specified grid.
//' @param X Double - Continuous value, to be identified on specified grid.
//' @param domain Vector of type double - Vector with two components containing
//' the bounds of the discretized interval
//' @param N_discr integer - discretization of the interval bounded by domain
//' @details This routine is implemented in C++. Note that a (continuous) value
//' between two grid
//' points condenses two its right grid point
//' @return Integer - The index position on the specified grid belonging to the continuous value x or -1
//' if input values are not adequate
//' @author P.v.W. Crommelin
//' @references S. Sprott, Chaos and Time-series analysis
//' @examples
//' // C++ C++ DEFINITION
//' double get_i(double x, Rcpp::NumericVector domain, int N_discr){
//'   if(N_discr>1.0 & domain[1]>domain[0] & x<=domain[1] & x>=domain[0]){
//'     double dx = (domain[1]-domain[0])/(N_discr-1.0);
//'     double minDiff = domain[1]-domain[0];
//'     int i = 0;
//'     double diff = abs(i*dx-x);
//'     while(diff <= minDiff){
//'       minDiff = diff;
//'       diff = abs(++i*dx - x);
//'     }
//'     return(i);
//'   }else{
//'     return(-1);
//'     std::printf("Input parameters not adequately chosen. See help(get_i) for more information");
//'   }
//' }
//'
//' /*** R
//' #R code example
//'
//' N_discr = 10 # Discretization of grid
//' lower = 0
//' upper = 1
//' domain = c(lower,upper)
//' dx = (domain[2]-domain[1])/(N_discr-1)
//' X = seq(0,1,dx) # grid (for comparism)
//' X2 = seq(0,1,dx/2) # test values to be indexed, note: step size is dx/2
//' (df = data.frame(TestValues = X2, CorrespondingIndexPositionsOnGrid = sapply(X2,function(x) get_i(x,domain,N_discr))))
//' */
//' @export
// [[Rcpp::export]]
double get_i(double x, Rcpp::NumericVector domain, int N_discr){
  if(N_discr>1.0 & domain[1]>domain[0] & x<=domain[1] & x>=domain[0]){
    double dx = (domain[1]-domain[0])/(N_discr-1.0);
    double minDiff = domain[1]-domain[0];
    int i = 0;
    double diff = abs(i*dx-x);
    while(diff <= minDiff){
      minDiff = diff;
      diff = abs(++i*dx - x);
    }
    return(i);
  }else{
    return(-1);
  }
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector vecMap_iter2(int N, double x0, Rcpp::NumericVector vecMap, Rcpp::NumericVector domain, int N_discr, bool skipFirst){
  Rcpp::NumericVector Series(N);
  if(skipFirst){
    Series[0] = vecMap[get_i(x0,domain,N_discr)];
  }else{
    Series[0] = x0;
  }
  for( int i=1; i<N ; i++){
    Series[i] = vecMap[get_i(Series[i-1],domain,N_discr)];
  }
  return(Series);
}

//' @export
// [[Rcpp::export]]
double vecMap(double x, Rcpp::NumericVector vec, Rcpp::NumericVector domain, int N_discr){
  double dx = (domain[1]-domain[0])/(N_discr-1.0);
  return(vec[get_i(x,domain,N_discr)]*dx);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector getIndexVector(Rcpp::NumericVector vec, Rcpp::NumericVector domain, int N_discr){
  Rcpp::NumericVector indices(N_discr);
  for(int j=0 ; j<vec.size() ; j++){
    indices[j] = get_i(vec[j],domain,N_discr);
  }
  return(indices);
}


//' @title general symmetric map
//' @description a generalization of a one dimensional map, which incorporates (among others) the lorenz-, logistic and cubeMap
//' @param x double - input value for one iteration
//' @param r double - control parameter
//' @param alpha double - exponent
//' @details This routine is implemented in C++
//' @return double - the result of one single iteration/
//' @author J.C. Lemm, P.v.W. Crommelin
//' @references S. Sprott, Chaos and Time-series analysis
//' @examples
//' //C++ C++ DEFINITION
//' double gsm_cpp(double x, double r, double alpha){
//'   return(r*(1-pow(2*abs(x-0.5),alpha)));
//' }
//' @export
// [[Rcpp::export]]
double gsm_cpp(double x, double r, double alpha){
  return(r*(1-pow(2*abs(x-0.5),alpha)));
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector indexGsm(double r, double alpha, int N_discr){
  Rcpp::NumericVector domain(2); //we are only interested for the general symmetric map mapping onto the interval [0,1]
                                //therefor hardcoded here
  domain[0] = 0.0;
  domain[1] = 1.0;
  double dx = (domain[1]-domain[0])/(N_discr-1.0);
  Rcpp::NumericVector codomain(N_discr); //evaluate gsm on grid points
  for(int i=0 ; i<N_discr ; i++){
    codomain[i] = gsm_cpp(i*dx,r,alpha);
  }
  //get indexposition on defined grid
  Rcpp::NumericVector indices = getIndexVector(codomain,domain,N_discr);
  return(indices);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector vecMap_iter(int N, double x0, Rcpp::NumericVector vec, Rcpp::NumericVector domain, int N_discr, bool skipFirst){
  Rcpp::NumericVector Series(N);
  if(skipFirst){
    Series[0] = vecMap(x0,vec,domain,N_discr);
  }else{
    Series[0] = x0;
  }
  for( int i=1 ; i<N ; i++){
    Series[i] = vecMap(Series[i-1],vec,domain,N_discr);
  }
  return(Series);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gsmAsVec_iter(int N, double x0, double r, double alpha, int N_discr, bool skipFirst){
  Rcpp::NumericVector domain(2); //we are only interested for the general symmetric map mapping onto the interval [0,1]
                                //therefor hardcoded here
  domain[0] = 0.0;
  domain[1] = 1.0;
  Rcpp::NumericVector gsmAsIndexVector = indexGsm(r,alpha,N_discr);
  return(vecMap_iter(N,x0,gsmAsIndexVector,domain,N_discr,skipFirst));
}



// private routine
double min_cpp(Rcpp::NumericVector vec){
  double m = vec[0];
  for (int i = 0 ; i<vec.size() ; i++){
    if(vec[i]<m){m = vec[i];}
  }
  return(m);
}

//' @title discretize the continuous general symmetric map
//' @description This function converts the continuous general symmetric map into a matrix
//' @param r double - control parameter
//' @param alpha double - exponent
//' @param N_discr integer - range of output matrix
//' @details This routine is implemented in C++.
//' @return matrix of type double - the corresponding matrix of range N
//' @author J.C. Lemm, P.v.W. Crommelin
//' @references S. Sprott, Chaos and Time-series analysis
//' @examples
//' Rcpp::NumericMatrix getMat(double r, double alpha, int N_discr){
//'   //grid is defined by N_discr values and therefore stepsize dx:
//'   double dx = 1.0/(N_discr-1.0);
//'   Rcpp::NumericVector codomain(N_discr);
//'   for(int i=0 ; i<N_discr ; i++){
//'     codomain[i] = gsm_cpp(i*dx,r,alpha);
//'   }
//'   //get indexposition on defined grid
//'   Rcpp::NumericVector index(N_discr);
//'   double diff(N_discr);
//'   for (int j=0 ; j<N_discr ; j++){
//'     double min_diff = 1.0; //1.0 is the maximum distance
//'     for (int i=0 ; i<N_discr ; i++){
//'       diff = abs(codomain[j] - i*dx);
//'       if(diff<=min_diff){
//'         min_diff = diff;
//'         index[j] = i;
//'         }
//'     }
//'   }
//'   Rcpp::NumericMatrix A(N_discr,N_discr);
//'   //consider "ROW MAJOR ORDER" (Wikipedia) to loop over arrays in c/c++
//'   //not implemented here
//'   for (int i=0 ; i<N_discr ; i++){
//'     for (int j=0 ; j<N_discr ; j++){
//'       if (j == index[i]){
//'         A(j,i) = 1.0;
//'       }else{
//'         A(j,i) = 0.0;
//'       }
//'     }
//'   }
//'   return(A);
//' }
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix getMat(double r, double alpha, int N_discr){
  //grid is defined by N_discr values and therefore stepsize dx:
  Rcpp::NumericVector index = indexGsm(r,alpha,N_discr);
  Rcpp::NumericMatrix A(N_discr,N_discr);
  //consider "ROW MAJOR ORDER" (Wikipedia) to loop over arrays in c/c++
  //not implemented here
  for (int i=0 ; i<N_discr ; i++){
    for (int j=0 ; j<N_discr ; j++){
      if (j == index[i]){
        A(j,i) = 1.0;
      }else{
        A(j,i) = 0.0;
      }
    }
  }
  return(A);
}

// private routine
Rcpp::NumericVector valToVec_cpp(double val, int N_discr, Rcpp::NumericVector domain){
  Rcpp::NumericVector gridDist(N_discr);
  double lower = domain[0];
  double upper = domain[1];
  double dx = (upper - lower)/(N_discr-1.0);
  for (int i=0; i<N_discr ; i++){
    gridDist[i] = abs(val-(i*dx+lower));
  }
  double min_gridDist = min_cpp(gridDist);
  int index = 0;
  for (int i=0; i<N_discr ; i++){
    if(gridDist[i]==min_gridDist){
      index = i;
    }
  }
  Rcpp::NumericVector result(N_discr);
  result[index] = 1;
  return(result);
}

// private routine
Rcpp::NumericVector dotProd(Rcpp::NumericMatrix A,Rcpp::NumericVector x){
  int n = A.nrow();
  int m = A.ncol();
  Rcpp::NumericVector b(n);
  if(m==x.size()){
    for (int i=0; i<n ; i++){
      for (int j=0; j<m ; j++){
        b[i] += A(i,j)*x[j];
      }
    }
  }
  return(b);
}

// private routine
// currently not in use
double innerProd(Rcpp::NumericVector a,Rcpp::NumericVector b){
  int n1 = a.size();
  int n2 = b.size();
  double result = 0;
  if(n1==n2){
    for (int i=0 ; i<n1 ; i++){
      result += a[i]*b[i];
    }
  }
  return(result);
}

// private routine
double vecToVal_cpp(Rcpp::NumericVector x,Rcpp::NumericVector domain){
  int n = x.size();
  double dx = (domain[1]-domain[0])/(n-1.0);
  double result = 0;
  for (int i=0 ; i<n ; i++){
    result += x[i]*(i*dx+domain[0]);
  }
  return(result);
}

//' @title discretize a real value
//' @description this routine discretizes a value corresponding to the desired method
//' @param val double - real value to be discretized
//' @param N_discr integer - desired discretization
//' @param domain vector of double - vector with two components containing the domain bounds
//' @param method integer - discretize with 1...round routine, or 2...custom method
//' @details This routine is implemented in C++
//' @return double - discretized value
//' @author J.C. Lemm, P. v.W. Crommelin
//' @examples
//' //C++ DEFINITION
//' double discretize_cpp(double val,  int N_discr, Rcpp::NumericVector domain, int method){
//'   double lower = domain[0];
//'   double upper = domain[1];
//'   if(method == 1){
//'     N_discr--;
//'     return(round(val*N_discr)/N_discr);
//'   }else{
//'     return(vecToVal_cpp(valToVec_cpp(val,N_discr,domain),domain));
//'   }
//' }
//' @export
// [[Rcpp::export]]
double discretize_cpp(double val,  int N_discr, Rcpp::NumericVector domain, int method){
  if(method == 1){
    N_discr--;
    return(round(val*N_discr)/N_discr);
  }else{
    return(vecToVal_cpp(valToVec_cpp(val,N_discr,domain),domain));
  }
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector iterFromMatrix_cpp(double x0, Rcpp::NumericMatrix mat, int N, Rcpp::NumericVector domain, bool skipFirst){
  int N_discr = mat.nrow();
  Rcpp::NumericVector x = valToVec_cpp(x0,N_discr,domain);
  Rcpp::NumericVector Series(N);
  //initialization
  if(skipFirst){
    x = dotProd(mat,x);
  }
  Series[0] = vecToVal_cpp(x,domain);
  for(int i=1 ; i<N ; i++){
    x = dotProd(mat,x);
    Series[i] = vecToVal_cpp(x,domain);
  }
  return(Series);
}

//' @title Time series creation with discretized output
//' @description Create time series produced by discrete model of the general symmetric map
//' @param N integer - Number of iterations
//' @param x0 double - starting value
//' @param r double - controll parameter
//' @param alpha double - exponent of general symmetric map
//' @param N_discr integer - controlls discretization of state space
//' @param skipFirst Boolean - If set to FALSE, the resulting time series contains the initial value x0
//' @details This routine is implemented in C++
//' @return vector of type double - the resulting time series
//' @author J.C. Lemm, P. v.W. Crommelin
//' @examples
//' //C++ DEFINITION
//' Rcpp::NumericVector Dgsm_iter_cpp(int N, double x0, double r, double alpha, int N_discr, bool skipFirst){
//'   Rcpp::NumericMatrix A = getMat(r,alpha,N_discr);
//'   //create vector from starting value
//'   Rcpp::NumericVector x = valToVec_cpp(x0,N_discr);
//'   Rcpp::NumericVector Series(N);
//'   //initialization
//'   if(skipFirst){
//'     x = dotProd(A,x);
//'   }
//'   Series[0] = vecToVal_cpp(x);
//'   for(int i=1 ; i<N ; i++){
//'     x = dotProd(A,x);
//'     Series[i] = vecToVal_cpp(x);
//'   }
//'   return(Series);
//' }
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector Dgsm_iter_cpp(int N, double x0, double r, double alpha, int N_discr, bool skipFirst){
  Rcpp::NumericMatrix A = getMat(r,alpha,N_discr);
  //create vector from starting value
  Rcpp::NumericVector domain(2);
  domain[0] = 0.0;
  domain[1] = 1.0;
  return(iterFromMatrix_cpp(x0,A,N,domain,skipFirst));
}

//' @title Time series creation
//' @description Create time series produced by general symmetric map
//' @param N integer - Number of iterations
//' @param x0 double - starting value
//' @param r double - controll parameter
//' @param alpha double - exponent of general symmetric map
//' @param N_discr integer - controlls discretization of state space (zero corresponds to continuous case)
//' @param skipFirst Boolean - If set to FALSE, the resulting time series contains the initial value x0
//' @param method integer - method of discretization with 1... round routine from base package (fast) or 2... matrix multiplication (slow)
//' @details This routine is implemented in C++
//' @return vector of type double - the resulting time series
//' @author J.C. Lemm, P. v.W. Crommelin
//' @examples
//' //C++ DEFINITION
//' Rcpp::NumericVector gsm_iter_cpp(int N, double x0, double r, double alpha, int N_discr,bool skipFirst, int method){
//'   if(method==1){
//'     if(N_discr>0){
//'       N_discr--; //Decrement because the rounding routine maps to N_discr + 1 values
//'     }
//'     Rcpp::NumericVector X(N);
//'     if(skipFirst){
//'       if(N_discr == 0){
//'         X[0] = gsm_cpp(x0,r,alpha);
//'       }else{
//'         X[0] = round(gsm_cpp(x0,r,alpha)*N_discr)/N_discr;
//'       }
//'     }else{
//'       X[0] = x0;
//'     }
//'     if(N>1){
//'       if(N_discr == 0){
//'         for(int i=1; i<N ; i++){
//'           X[i] = gsm_cpp(X[i-1],r,alpha);
//'         }
//'       }else{
//'         for(int i=1; i<N ; i++){
//'           X[i] = round(gsm_cpp(X[i-1],r,alpha)*N_discr)/N_discr;
//'         }
//'       }
//'     }
//'     return(X);
//'   }else{
//'     return(Dgsm_iter_cpp(N,x0,r,alpha,N_discr,skipFirst));
//'   }
//' }
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gsm_iter_cpp(int N, double x0, double r, double alpha, int N_discr,bool skipFirst, int method){
  if(method==1){
    if(N_discr>0){
      N_discr--; //Decrement because the rounding routine maps to N_discr + 1 values
    }
    Rcpp::NumericVector X(N);
    if(skipFirst){
      if(N_discr == 0){
        X[0] = gsm_cpp(x0,r,alpha);
      }else{
        X[0] = round(gsm_cpp(x0,r,alpha)*N_discr)/N_discr;
      }
    }else{
      X[0] = x0;
    }
    if(N>1){
      if(N_discr == 0){
        for(int i=1; i<N ; i++){
          X[i] = gsm_cpp(X[i-1],r,alpha);
        }
      }else{
        for(int i=1; i<N ; i++){
          X[i] = round(gsm_cpp(X[i-1],r,alpha)*N_discr)/N_discr;
        }
      }
    }
    return(X);
  }else{
    return(Dgsm_iter_cpp(N,x0,r,alpha,N_discr,skipFirst));
  }
}

//' @title Gaussian likelihood for given time series data
//' @description Lik_gsm_cpp produces a time Series generated by the general symmetric map (gsm)
//' and computes the its gaussian likelihood given the supplied data. The standart deviation sigma is
//' equal for every datapoint in this implementation. Eventually supply a single scalar value of type double.
//' @param alpha Double - Exponent of the gsm.
//' @param r Double - Control parameter of the gsm.
//' @param x0 Double - Initial value of the time series.
//' @param Y Vector of type double - The given data.
//' @param sigma Double - The standard deviation of the gaussian likelihood.
//' @param N_discr Integer - Discretization of the state space (N_discr = 0 -> continuous case).
//' @param method Integer - Discretization with 1... R's core round routine (fast) or 2...A Matrix multiplication (slower)
//' @details This routine is implemented in C++.
//' @author J.C. Lemm, P. v.W. Crommelin
//' @references S. Sprott, Chaos and Time-series analysis
//' @examples
//' r = 0.99 # True but hidden control parameter of the general symmetric map (to be reconstructed)
//' x0 = 0.2 # starting value
//' alpha = 2.0
//' N = 20
//' SIGMA = 0.1
//' N_discr = 0 # coded value for continuous time series
//'   Series = rep(x0,N)
//'   Series = myBayes::gsm_iter_cpp(alpha = alpha,
//'                                  r = r,
//'                                  x0 = x0,
//'                                  N = N,
//'                                  skipFirst = TRUE,
//'                                  N_discr = N_discr,
//'                                  method = 1)
//'   Y = Series + rnorm(n = N,
//'                      mean = 0,
//'                      sd = SIGMA) # add noise (gauss distributed)
//'
//'   vec = seq(from = 0.95,
//'             to = 1.0,
//'             by = 0.000001)
//'
//'   Lik = sapply(X = vec,
//'                FUN = function(val) myBayes::Lik_gsm_cpp(alpha = alpha,
//'                               r = val,
//'                               x0 = x0,
//'                               Y = Y,
//'                               sigma = SIGMA,
//'                               N_discr = N_discr,
//'                               method = 1))
//'   par(mfrow = c(2,1))
//'   plot(x = 1:N,
//'        y = Y,
//'        main = "Given data",
//'        xlab = "iteration",
//'        ylab = "Value of timeSeries",
//'        col = "red",
//'        type = "l")
//'   points(x = 1:N,
//'          y = Y,
//'          col = "red",
//'          cex = 1.5,
//'          pch = 16)
//'   grid()
//'   plot(x = vec,
//'        y = Lik,
//'        main = "Gaussian likelihood as a function of control parameter r given the data",
//'        xlab = "control parameter r",
//'        ylab = "Likelihood",
//'        cex = 0.5,
//'        pch = 16)
//'   lines(x = rep(r,2),
//'         y = c(0,max(Lik)),
//'         col = "red",
//'         lwd = 3,
//'         lty = 1)
//'   grid()
//'   legend("topleft",
//'          leg = paste0("Real value for r = ",r),
//'          col = "red",
//'          lty = 1,
//'          lwd = 3)
//'
//'
//' #C++ DEFINITION:
//' #double Lik_gsm_cpp(double alpha, double r, double x0, NumericVector Y, double sigma, int N_discr){
//' #  int n = Y.size(); //Because equally the generated skips first datapoint
//' #  bool skipFirst = true;
//' #  Rcpp::NumericVector X = gsm_iter_cpp(n,x0,r,alpha,N_discr,skipFirst);
//' #  double sum = 0;
//' #  for(int i = 0; i<n ; i++){
//' #    sum += pow(Y[i]-X[i],2.0);
//' #  }
//' #  double L = pow(2.0*PI*pow(sigma,2),-0.5*n)*exp(-0.5*sum/pow(sigma,2.0));
//' #  return(L);
//' #}
//' @export
// [[Rcpp::export]]
double Lik_gsm_cpp(double alpha, double r, double x0, Rcpp::NumericVector Y, double sigma, int N_discr, int method){
  int n = Y.size();
  Rcpp::NumericVector X(n);
  if(N_discr==0){
    X = gsm_iter_cpp(n,x0,r,alpha,N_discr,true,method);
  }else{
    if(method==1){
      X = gsm_iter_cpp(n,x0,r,alpha,N_discr,true,method);
    }else{
      X = Dgsm_iter_cpp(n,x0,r,alpha,N_discr,true);
    }
  }
  double sum = 0;
  for(int i = 0; i<n ; i++){
    sum += pow(Y[i]-X[i],2.0);
  }
  //double L = pow(2.0*PI*pow(sigma,2),-0.5*n)*exp(-0.5*w*sum/pow(sigma,2.0));
  double L = exp(-0.5*sum/pow(sigma,2.0));
  return(L);
}

//' @title Gaussian likelihood for given data
//' @description Generates a time Series generated by the DISCRETIZED version of the general symmetric map (gsm) and computes the likelihood given the data
//' @param alpha double - exponent of gsm
//' @param r double - control parameter of gsm
//' @param x0 double - initial value of the time series
//' @param Y vector of type double - the given data
//' @param sigma double - the standard deviation of the gaussian likelihood
//' @param N_discr integer - discretization of the state space
//' @details This routine is implemented in C++
//' @author J.C. Lemm, P. v.W. Crommelin
//' @references S. Sprott, Chaos and Time-series analysis
//' @examples
//' //C++ DEFINITION:
//' double Lik_Dgsm_cpp(double alpha, double r, double x0, Rcpp::NumericVector Y, double sigma, int N_discr){
//'   int n = Y.size();
//'   bool skipFirst = true;
//'   Rcpp::NumericVector X = Dgsm_iter_cpp(n,x0,r,alpha,N_discr,skipFirst);
//'   double sum = 0;
//'   for(int i = 0; i<n ; i++){
//'     sum += pow(Y[i]-X[i],2.0);
//'   }
//'   double L = pow(2.0*PI*pow(sigma,2),-0.5*n)*exp(-0.5*sum/pow(sigma,2.0));
//'   return(L);
//' }
//' @export
// [[Rcpp::export]]
double Lik_Dgsm_cpp(double alpha, double r, double x0, Rcpp::NumericVector Y, double sigma, int N_discr){
  int n = Y.size();
  Rcpp::NumericVector X = Dgsm_iter_cpp(n,x0,r,alpha,N_discr,true);
  double sum = 0;
  for(int i = 0; i<n ; i++){
    sum += pow(Y[i]-X[i],2.0);
  }
  double L = pow(2.0*PI*pow(sigma,2),-0.5*n)*exp(-0.5*sum/pow(sigma,2.0));
  return(L);
}

