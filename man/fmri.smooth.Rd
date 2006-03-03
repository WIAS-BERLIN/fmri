\name{fmri.smooth}
\alias{fmri.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothing Statistical Parametric Maps }
\description{
  Perform the adaptive weights smoothing procedure 
}
\usage{
fmri.smooth(spm, hmax = 4, type = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ object of class \code{fmri.spm} }
  \item{hmax}{ maximum bandwidth to smooth }
  \item{type}{ if "old" switch to older version of AWS. for advanced use
    only! }
  \item{\dots}{ Other arguments passed to vaws3D }
}
\details{
  This function is a wrapper for the vaws3D function. It first
  eliminates zero variances by assigning 1e20 (thus will be excluded
  from smoothing!)

  \code{hmax} is the maximum bandwidth of the iteration. \code{qlambda} is the
  adaptation parameter. Choose 1 (and hinit > hmax) for non adaptive
  smoothing here use default  value for adaptive smoothing. \code{lkern} can be used for specifying the
  localization kernel. For comparison with non adaptive methods use
  "Gaussian" (hmax given in FWHM), for better adaptation use "Triangle"
  (hmax given in voxel). See \code{\link{vaws3D}} for more details).

}
\value{
  object of class \code{fmrispm}
  \item{cbeta}{ smoothed parameter estimate}
  \item{var}{ variance of the parameter}
  \item{hmax}{ achieved maximum bandwidth }
  \item{rxyz}{ smoothness in resel space. all directions }
  \item{scorr}{ spatial correlation of original data }
  \item{weights}{ ratio of voxel dimensions }
  \item{smooth}{ indicates that SPM is smoothed. "TRUE" }
  \item{vwghts}{ ratio of estimated variances for the stimululi given by
  \code{vvector}}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{vaws3D}} }
\examples{
\dontrun{fmri.smooth(spm, hmax = 4, lkern = "Gaussian", ...)}
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
