\name{get3Dh.gauss}
\alias{get3Dh.gauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothness }
\description{
  Get corresponding bandwidth for local variance reduction.
}
\usage{
get3Dh.gauss(vred, h0, vwghts, step = 1.002, interv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vred}{ variance reduction }
  \item{h0}{ vector of bandwidth with  }
  \item{vwghts}{ vector of weights }
  \item{step}{ step size for iteration }
  \item{interv}{ discretization parameter for gauss kernel }
}
\details{
  Calculate the bandwidth corresponding to the amount of variance
  reduction achieved locally.
}
\value{
  bandwidth
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{vaws3D}} }
\examples{
\dontrun{get3Dh.gauss(.01, h0=c(.5,.7,.3), vwghts=c(1,1,1), step = 1.002, interv = 1)}
}
\keyword{internal }% at least one, from doc/KEYWORDS
