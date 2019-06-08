#### general symmetric map (C extension) ####
#sig = c(N="integer",x0="numeric",r="numeric",alpha="numeric",skipFirst="integer")
#body = "
##include <math.h>
##include <Rinternals.h>
#
#
#// convert every SEXP object into a c datatype
#double r_ = asReal(r);
#double a_ = asReal(alpha);
#double x0_ = asReal(x0);
#int N_ = asInteger(N);
#int skipFirst_ = asInteger(skipFirst);
#
#
#// Allocate an output vector
#SEXP X = PROTECT(allocVector(REALSXP,N_));
#
#// create an Array pointing to that output vector
#double *X_;
#X_ = REAL(X);
#
#
#//Manipulate the array
#if (skipFirst_)  X_[0] = r_*(1.0 - pow(2.0 * fabs(x0_-0.5),a_));
#if (!skipFirst_) X_[0] = x0_;
#
#int i;
#for ( i=0 ; i<(N_-1) ; i++ ){
#    X_[i+1] = r_*(1.0 - pow(2.0 * fabs(X_[i]-0.5),a_));
#}
#
#UNPROTECT(1);
#//return that array
#return X;
#"
#
#gensymMap_iter_c =  inline::cfunction(sig = sig,
#                              body = body,
#                              verbose = FALSE,
#                              convention = ".Call",
#                              language = "C")
#setCMethod(f = "gensymMap_iter_c",
#           sig = sig,
#           body = body,
#           verbose = FALSE)

#test = inline::cfunction(sig = c(x="integer"),
#                         body = 'return(x);')
