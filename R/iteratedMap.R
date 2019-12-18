#' @title Time series creation by given map.
#' @param x0 double - initial condition.
#' @param n integer - length of resulting time series.
#' @param fun function of one parameter - the corresponding map.
#' @param skipfirst logical - wether the initial condition is to be removed.
#' @details This routine takes a map and applies it iteratively.
#' @return vector of type double and length n - the resulting time series.
#' @author Philipp van Wickevoort Crommelin
#' @examples
#'N = 12
#'x = seq(from = 0,
#'        to = 1,
#'        by = 1/(N-1))
#'y = stats::runif(N)
#'fun = stats::approxfun(x = x,
#'                       y = y,
#'                       method = "linear")
#'x0 = 0.2
#'par(mfrow = c(2,1))
#'plot(x = x,
#'     y = y,
#'     main = "map",
#'     type = "l",
#'     col = "lightblue",
#'     lwd = 7,
#'     ylab = "")
#'N_series = 30
#'tmser = myBayes::iteratedMap(n = N_series,
#'                    x = x0,
#'                    fun = fun,
#'                    skipFirst = FALSE)
#'plot(NULL,
#'     xlim = c(1,N_series),
#'     ylim = c(1,0),
#'     main = "Time series",
#'     xlab = "iteration",
#'     ylab = "")
#'lines(x = 1:N_series,
#'      y = tmser,
#'      col = "purple",
#'      lwd = 5,
#'      ylab = "")
#'points(x = 1:N_series,
#'       y = tmser,
#'       col = "orange",
#'       cex = 2,
#'       pch = 16)
#' @export
iteratedMap = function(n,x0,fun,skipFirst = TRUE){
  X = rep(0,n)
  if(skipFirst){
    X[1] = fun(x0)
  }else{
    X[1] = x0
  }
  if(n>=2){
    for(i in 2:n){
      X[i] = fun(X[i-1])
    }
  }
  return(X)
}



