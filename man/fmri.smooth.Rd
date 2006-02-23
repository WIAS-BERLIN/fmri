\name{fmri.smooth}
\alias{fmri.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothing Statistical Parametric Maps }
\description{
  Perform the adaptive weights smoothing procedure 
}
\usage{
fmri.smooth(spm, hmax = 4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ object of class \code{fmri.spm} }
  \item{hmax}{ maximum bandwidth to smooth }
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
  object of class \code{fmri.spm}
  \item{cbeta}{ smoothed parameter estimate}
  \item{var}{ variance of the parameter}
  \item{hmax}{ achieved maximum bandwidth }
  \item{rxyz}{ smoothness in resel space. all directions }
  \item{scorr}{ spatial correlation of original data }
  \item{weights}{ ratio of voxel dimensions }
  \item{smooth}{ indicates that SPM is smoothed. "TRUE" }
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{vaws3D}}, \code{\link{correlation}} }
\examples{
\dontrun{fmri.smooth(spm, hmax = 4, lkern = "Gaussian", ...)}
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
