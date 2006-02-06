\name{Varcor.gauss}
\alias{Varcor.gauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothness }
\description{
  Calculates a correction for the variance estimate
}
\usage{
Varcor.gauss(h, interv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{ bandwidth }
  \item{interv}{ discretization parameter for gauss kernel }
}
\details{
  Calculates a correction for the variance estimate obtained by
  \eqn{(IQRdiff(y)/1.908)^2}{(IQRdiff(y)/1.908)^2} in case of colored noise that was produced by
  smoothing with gaussian kernel and bandwidth \code{h}.

  \code{interv} allows for further discretization of the Gaussian
  Kernel, result depends on
  \code{interv} for small bandwidths. \code{interv}=1  is correct for kernel smoothing, 
  \code{interv}>>1 should be used to handle intrinsic correlation (smoothing preceeding 
  discretisation into voxel) 
}
\value{
  Variance estimate correction factor.
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{Spatialvar.gauss}} }
\examples{
  # Example 1
  Varcor.gauss(3, interv = 1)
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
