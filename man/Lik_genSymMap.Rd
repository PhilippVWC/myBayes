\name{Lik_genSymMap}
\alias{Lik_genSymMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Gaussian Likelihood for general symmetric map
}
\description{
  The gaussian likelihood is a metric which describes how good a hypothesis matches a given set of observational data.
}
\usage{
  Lik_genSymMap(alpha, r, x0, Y, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
  exponent
}
  \item{r}{
  control parameter
}
  \item{x0}{
  starting value
}
  \item{Y}{
  given set of Data
}
  \item{sigma}{
  the standard deviation of the underlying gaussian likelihood
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a scalar value for the likelihood
}
\references{
S. Sprott, Chaos and Time-series analysis}
\author{
J.C. Lemm, P. van Wickevoort Crommelin
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## The function is currently defined as
function (alpha, r, x0, Y, sigma)
{
    n = length(Y)
    X = gensymMap_iter_c(N = n, x0 = x0, r = r, alpha = alpha)
    L = LLik(X = X, Y = Y, sigma = sigma)
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ likelihood }% use one of  RShowDoc("KEYWORDS")
\keyword{ gauss }% __ONLY ONE__ keyword per line
