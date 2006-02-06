\name{Spatialvar.gauss}
\alias{Spatialvar.gauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothness }
\description{
  Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h
}
\usage{
Spatialvar.gauss(h, h0, d, interv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{ bandwidth }
  \item{h0}{ bandwidth }
  \item{d}{ dimension of the data }
  \item{interv}{ discretization parameter for gauss kernel }
}
\details{
  Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h in 
  case of colored noise that was produced by smoothing with Gaussian kernel and bandwidth h0
  
  Spatialvar.gauss(h,h0,d)/Spatialvar.gauss(h,1e-5,d) gives the 
  a factor for lambda to be used with bandwidth h 
}
\value{factor of variance reduction
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{Varcor.gauss}} }
\examples{
  # Example 1
  Spatialvar.gauss(4,1e-5,3)
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
