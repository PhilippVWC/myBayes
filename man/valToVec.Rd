\name{valToVec}
\alias{valToVec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Translate a scalar value to a probability vector
}
\description{
This routine translates a scalar value into a vector of length N containing N-1 zeros and one 1.0.
The index position of the 1.0 within the returning vector depends on the position of val with respect to
the bounds x_lower and x_upper.
}
\usage{
valToVec(val, N = 10, x_lower = 0, x_upper = 1)
}
\arguments{
  \item{val}{
  Value to be translated
}
  \item{N}{
  Dimension of resulting vector
}
  \item{x_lower}{
  Lower bound
}
  \item{x_upper}{
  Upper bound
}
}
\details{
Always provide N>1. If you convert the resulting vector back with <<vecToVal>> (see Package "myBayes"). The resulting value deviates from
the original by a maximum difference of +- 1/N.
}
\value{
A probability vector containing N-1 zeros and one 1.0.
}
\references{
%% ~put references to the literature/web site here ~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (val, N = 10, x_lower = 0, x_upper = 1)
{
    if (val > x_upper | val < x_lower)
        print("warning: val is not within bounds of x_lower and x_upper")
    if (N > 1) {
        dx = (x_upper - x_lower)/(N - 1)
        ref = seq(from = x_lower, to = x_upper, by = dx)
        diff = abs(ref - rep(val, N))
        minInd = which(diff == min(diff))
        result = rep(0, N)
        result[minInd] = 1
        return(result)
    }
    else {
        print(paste0("N needs to be greater than 1. N is currently set to ",
            N))
        return(NA)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
