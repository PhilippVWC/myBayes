# r = 0.99 # True but hidden control parameter of the general symmetric map (to be reconstructed)
# x0 = 0.2 # starting value
# alpha = 2.0
# N = 20
# SIGMA = 0.1
# N_discr = 0 # coded value for continuous time series
# Series = rep(x0,N)
# Series = myBayes::gsm_iter_cpp(alpha = alpha,
#                                r = r,
#                                x0 = x0,
#                                N = N,
#                                skipFirst = TRUE,
#                                N_discr = N_discr,
#                                method = 1)
# Y = Series + rnorm(n = N,
#                    mean = 0,
#                    sd = SIGMA) # add noise (gauss distributed)
#
# vec = seq(from = 0.95,
#           to = 1.0,
#           by = 0.000001)
#
# Lik = sapply(X = vec,
#              FUN = function(val) myBayes::Lik_gsm_cpp(alpha = alpha,
#                                                       r = val,
#                                                       x0 = x0,
#                                                       Y = Y,
#                                                       sigma = SIGMA,
#                                                       N_discr = N_discr,
#                                                       method = 1))
# par(mfrow = c(2,1))
# plot(x = 1:N,
#      y = Y,
#      main = "Given data",
#      xlab = "iteration",
#      ylab = "Value of timeSeries",
#      col = "red",
#      type = "l")
# points(x = 1:N,
#        y = Y,
#        col = "red",
#        cex = 1.5,
#        pch = 16)
# grid()
# plot(x = vec,
#      y = Lik,
#      main = "Gaussian likelihood as a function of control parameter r given the data",
#      xlab = "control parameter r",
#      ylab = "Likelihood",
#      cex = 0.5,
#      pch = 16)
# lines(x = rep(r,2),
#       y = c(0,max(Lik)),
#       col = "red",
#       lwd = 3,
#       lty = 1)
# grid()
# legend("topleft",
#        leg = paste0("Real value for r = ",r),
#        col = "red",
#        lty = 1,
#        lwd = 3)
#
