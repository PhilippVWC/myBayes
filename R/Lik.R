#' @title Gaussian likelihood for given time series data
#' @description Lik computes the Gaussian likelihood for a the supplied data Y given the model vector X.
#' The standart deviation sigma is
#' equal for every datapoint in this implementation. Hence, please supply a single scalar value of type double.
#' @param Y Vector of type double - The given data.
#' @param X Vector of type double - The given model vector.
#' @param sigma Double - The standard deviation of the gaussian likelihood.
#' @return Double - The resulting likelihood.
#' @author Philipp van Wickevoot Crommelin
#' @export
Lik = function(X,Y,sigma){
  return(exp(-sum(X-Y)^2/(2*sigma^2)))
}
