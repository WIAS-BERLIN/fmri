\name{fmri.design}
\alias{fmri.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model for FMRI Data }
\description{
  Returns a designmatrix for a linear model with given stimuli and
  possible polynomial drift terms.
}
\usage{
fmri.design(hrf, order = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hrf}{ matrix containing expected BOLD repsonse(s) for the linear model }
  \item{order}{ order of the polynomial drift terms }
}
\details{
  The stimuli given in \code{hrf} are used as first columns in
  the designmatrix.

  The order of the polynomial drift terms is given
  by \code{order}, which defaults to 1.

  The polynomials are defined orthogonal to the stimuli given in
  \code{hrf}. 
}
\value{
  designmatrix of the linear model
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{fmri.stimulus}}, \code{\link{fmri.lm}} }
\examples{
  # Example 1
  hrf <- fmri.stimulus(107, c(18, 48, 78), 15, 2)
  z <- fmri.design(hrf,2)
  plot(hrf,type="l")
}
\keyword{ regression }% at least one, from doc/KEYWORDS
\keyword{ design }% __ONLY ONE__ keyword per line
