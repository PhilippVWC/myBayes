\name{max2d}
\alias{max2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the maxima of a 2D matrix.
}
\description{
Function for 2D matrices returning extremal values and its corresponding indices.
}
\usage{
max2d = function(mat,epsilon,maximum=TRUE,maxRows=NA)
}
%- maybe also 'usage' for other objects documented here.
  \arguments{
  \item{epsilon}{
Relative tolerance level in [0,1], controlling the amount of extremal values desired.
Zero: consider only maximum values (possibly degenerated)
One: consider every value
}
\item{mat}{
Matrix at hand.
}
\item{maximum}{
Boolean value indicating, wether to search for a minimum or a maximum.
}
\item{maxRows}{
Upper bound for number of results.
}
}
\details{
If maxRows is provided, the returning data frame is sorted by MaxVal.
}
\value{
Data frame with three columns:
\item{MaxVal}{Extremal value.}
\item{RowIndex}{Corresponding x index of the extremum within the matrix.}
\item{ColIndex}{Corresponding y index of the extremum within the matrix.}
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
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
