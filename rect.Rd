\name{rect}
\alias{rect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rectangle function
}
\description{
produces a discontinuous function which is zero everywhere except within [a,b]
}
\usage{
rect(x, a, b, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  variable
}
  \item{a}{
  left bound
}
  \item{b}{
  right bound
}
  \item{z}{
  Height of the rectangle
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  function evaluation of x
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
function (x, a, b, z)
{
    result = z * (Heaviside(x = x, x0 = a) - Heaviside(x = x,
        x0 = b))
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
