\name{discretizeMap}
\alias{discretizeMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discretize a continuous function by a specified level
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
discretizeMap(N, map, xlim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{map}{
%%     ~~Describe \code{map} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (N, map, xlim)
{
    x_lower = xlim[1]
    x_upper = xlim[2]
    dx = (x_upper - x_lower)/(N - 1)
    domain = seq(from = x_lower, to = x_upper, by = dx)
    codomain = sapply(X = domain, FUN = function(x) map(x))
    helperMatrix = matrix(codomain, ncol = 1) \%*\% rep(1, N)
    indices = apply(X = helperMatrix, MARGIN = 1, FUN = function(x) {
        dist = abs(x - domain)
        index = which(dist == min(dist))
        return(index)
    })
    A = matrix(rep(0, N * N), ncol = N, nrow = N)
    for (i in 1:N) {
        A[indices[i], i] = 1
    }
    return(A)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
