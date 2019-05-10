\name{MH}
\alias{MH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampling from a custom distribution
}
\description{
Metropolis Hastings Algoritm. Create a random variable, which is distributed according to a custom probability distribution provided.
}
\usage{
MH(T = 1, x0, trKern = dnorm, sampler, N = 1000, bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
  "Temperature". This value determines the jumping distance for new candidate values. Choose within [0,1]
}
  \item{x0}{
  Starting Value. Choose within bounds.
}
  \item{trKern}{
  Transition Kernel. Function of at least to parameters. 1st param: new candidate, 2nd param: last position, corresponds to maximum of trKern. Default to normal distribution (symmetric).
  If you choose an asymmetric continuous probability distribution you impose a preferred sampling direction to the left or to the right.
}
  \item{sampler}{
  Custom probability distribution to draw samples from.
}
  \item{N}{
  Number of samples.
}
  \item{bounds}{
  Vector with two components constraining the domain. 1st comp: left bound, 2nd comp: right bound.
}
}
\details{
An asymmetric transition Kernel like a logitnormal distribution, beta distribution or a triangular distribution centered on the left, causes the
algorithm to prefer new candidates to the right of the current value and vice versa.
}
\value{
A vector of size N containing the samples.
}
\references{
Rubinstein, Simulation and the Monte Carlo method.
}
\author{
Philipp van Wickevoort Crommelin
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# ---------------------------------
# needs package "logitnorm"

# example without transition kernel
x = seq(0,1,0.01)
x0 = 0.4
bounds = c(0,1)
N = 50
mu = -1
sampler = function(x) dlogitnorm(x = x,
                                 mu = mu,
                                 sigma = 1)
MC = MH(x0 = x0,
        bounds = bounds,
        sampler = sampler)
par(mfrow = c(2,1))
plot(x,sampler(x))
hist(MC)




# Example with tranistion Kernel
projectUnityToInf = function(x) return( tan((x*pi)-pi/2) )
trKern = function(y,x) max(1e-14,dlogitnorm(x = y,
                                  mu = projectUnityToInf(x)))
#y = sapply(X = x,
#           FUN = function(x_i) trKern(y = x_i,
#                                      x = 0.8))
#par(mfrow=c(1,1))
#plot(x = x,
#     y = y,
#     main = paste0("Example transition Kernel"))

MC = MH(x0 = x0,
        bounds = bounds,
        sampler = sampler,
        trKern = trKern)

par(mfrow = c(2,1))
plot(x = x,
     y = sampler(x),
     ylab = "sampler")
hist(MC)
#---------------------------------
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
