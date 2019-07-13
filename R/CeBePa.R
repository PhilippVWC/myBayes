#' @title Parameter creation for Beta distribution.
#' @description This function transforms the mean value and the standart deviation of a univariate normal distribution
#' into two values alpha and beta, which can be used for a similar shaped distribution function, that is defined
#' on a bounded domain like the Beta distribution.
#' @param mean Double - A value in the interval [0,1].
#' @param sigmaRel Double - A value describing the width of the resulting maximum value.
#' @return A vector with two components of type double - The parameters alpha and beta.
#' @author J.C. Lemm, P.v.W. Crommelin
#' @examples
#' mean = 0.5
#' sigmaRel = 0.1
#'
#' x = seq(from = -10,
#'         to = 10,
#'         by = 0.01)
#' plot(x = x,
#'      y = dnorm(x = x,
#'                mean = mean,
#'                sd = sigmaRel),
#'      xlim = c(-0.1,1.1),
#'      ylim = c(0,10),
#'      type = "l")
#' x2 = seq(from = 0,
#'          to = 1,
#'          by = 0.001)
#' params = CeBePa(mean = mean,
#'                 sigmaRel = sigmaRel)
#' lines(x = x2,
#'       y = dbeta(x = x2,
#'                 shape1 = params[1],
#'                 shape2 = params[2]),
#'       col = "red",
#'       lwd = 3)
#' legend("topright",
#'        leg = c("Normal distribution","Beta distribution"),
#'        col = c("black","red"),
#'        lwd = c(1,3))
#' @export
CeBePa = function(mean,sigmaRel){
  N = 1/sigmaRel^2 - 1.0
  alpha = max(1, mean*N)
  beta = max(1, (1-mean)*N)
  return(c(alpha,beta))
}
