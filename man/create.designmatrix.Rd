\name{create.designmatrix}
\alias{create.designmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model for FMRI Data }
\description{
  Returns a designmatrix for a linear model with given stimulus and
  possible polynomial drift terms.
}
\usage{
create.designmatrix(hrf, order = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hrf}{ matrix containing signal functions for the linear model }
  \item{order}{ order of the polynomial drift terms }
}
\details{
  The stimuli given in \code{hrf} are used as first columns in
  the designmatrix.

  The order of the polynomial drift terms is given
  by\code{order}, which defaults to 1.

  The polynomials are defined orthogonal to the stimuli given in
  \code{hrf}. 
}
\value{
  designmatrix of the linear model
}
\author{ Karsten Tabelow (tabelow@wias-berlin.de) }
\seealso{ \code{\link{create.stimulus}}, \code{\link{calculate.lm}} }
\examples{
  # Example 1
  hrf <- create.stimulus(107, c(18, 48, 78), 15, 2)
  z <- create.designmatrix(hrf,2)
  dim(hrf)
}
\keyword{ regression }% at least one, from doc/KEYWORDS
\keyword{ design }% __ONLY ONE__ keyword per line
