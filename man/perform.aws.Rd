\name{perform.aws}
\alias{perform.aws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothing Statistical Parametric Maps }
\description{
  Perform the adaptive weights smoothing procedure 
}
\usage{
perform.aws(beta, variance, hmax = 4, hinit = 1, weights = c(1, 1, 1), qlambda = 1, lkern = "Gaussian", scorr = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{ 3 or 4 dimensional parameter map to smooth }
  \item{variance}{ variance estimate of the parameter map. }
  \item{hmax}{ maximum bandwidth to smooth }
  \item{hinit}{ initial bandwidth of the smoothing iteration. }
  \item{weights}{ weights vector for smoothing describing slice
    distances in x-, y-, and z-direction. }
  \item{qlambda}{ parameter qlambda for AWS procedure }
  \item{lkern}{ localization kernel to be used }
  \item{scorr}{ vector with estimated correlation in x-, y-, and z-direction. }
  \item{\dots}{ Other arguments passed to vaws3D }
}
\details{
  This function is a wrapper for the vaws3D function. It first
  eliminates zero variances by assigning variances below 25-percent
  quantile to 25-percent quantile.

  \code{hmax} is the maximum bandwidth of the iteration. \code{hinit} is
  the minimum bandwidth of the iteration. \code{qlambda} is the
  adaptation parameter. Choose 1 (and hinit > hmax) for non adaptive
  smoothing here. \code{lkern} can be used for specifying the
  localization kernel. For comparison with non adaptive methods use
  "Gaussian", for better adaptation use "Triangle" (see
  \code{\link{vaws3D}} for more details).

  \code{scorr} is the estimated correlation in x-, y-, and z-direction.
}
\value{
  \item{hat}{ smoothed parameter estimate}
  \item{var}{ variance of the parameter}
  \item{ni}{ sum of weights}
  \item{hmax}{ achieved maximum bandwidth }
  \item{lseq}{ The sequence specified in \code{lseq} of vaws3D }
  \item{y}{ The data. }
  \item{mae}{ MAE }
  \item{vred}{ achieved variance reduction}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{vaws3D}}, \code{\link{correlation}} }
\examples{
\dontrun{perform.aws(beta, variance, hmax = 4, hinit = 1, weights = c(1, 1, 1), qlambda = 1, lkern = "Gaussian", scorr = 0, ...)}
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
