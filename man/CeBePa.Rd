\name{CeBePa}
\alias{CeBePa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Centralize Beta Distribution Parameters alpha and beta
}
\description{
This function takes a relative standard deviation sigmaRel and mean and produces a beta distribution centered around that mean value.
}
\usage{
CeBePa(mean,sigmaRel)
}
\arguments{
  \item{mean}{
  Mean value, which will correspond to the most probable value of the resulting beta distribution
}
  \item{sigmaRel}{
  Relative "standard deviation". This value determines the steepness of the resulting beta distribution.
  Low values correspond to a sharp peak with high steepness and vice versa.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a vector with two components.
\item{1st component }{alpha}
\item{2nd Component }{beta}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jörg C. Lemm, Philipp van Wickevoort Crommelin
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
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
