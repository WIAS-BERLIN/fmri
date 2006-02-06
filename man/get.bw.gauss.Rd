\name{get.bw.gauss}
\alias{get.bw.gauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothness }
\description{
  Get the   bandwidth for lkern corresponding to a given correlation
}
\usage{
get.bw.gauss(corr, step = 1.001, interv = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corr}{ correlation }
  \item{step}{ step size for iteration }
  \item{interv}{ discretization parameter for gauss kernel }
}
\details{
  get the   bandwidth for lkern corresponding to a given correlation
  keep it simple result does not depend on d

  \code{interv} allows for further discretization of the Gaussian
  Kernel, result depends on
  \code{interv} for small bandwidths. \code{interv}=1  is correct for kernel smoothing, 
  \code{interv}>>1 should be used to handle intrinsic correlation (smoothing preceeding 
  discretisation into voxel) 
}
\value{
  bandwidth
}
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de} }
\examples{
  # Example 1
  get.bw.gauss(0.8)
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
