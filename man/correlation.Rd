\name{correlation}
\alias{correlation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothness }
\description{
  Estimate correlation
}
\usage{
correlation(res, mask = array(1, dim = dim(res)[1:3]))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{ array of residuals }
  \item{mask}{ logical. mask which area should be included in estimation }
}
\details{
  Estimate the correlation in x-, y-, z- direction restricted to
  \code{mask}.
}
\value{
  Vector of correlation in x-, y-, z- direction.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
  # Example 1
  correlation(array(rnorm(32*32*32*20),c(32,32,32,20)))
}
\keyword{ misc }% at least one, from doc/KEYWORDS
