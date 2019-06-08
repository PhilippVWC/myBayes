## (see Abramowitz and Stegun 29.2.29)
## 'complementary error function'
#' @export
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
