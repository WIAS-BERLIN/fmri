\name{fmri.plot}
\alias{fmri.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  Plot the fMRI detection result.
}
\usage{
fmri.plot(detect, anatomic = array(0, dim = dim(detect$detect)), pos = c(-1, -1, -1), type = "slice", device = "X11", file = "plot.png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{detect}{ object of class \code{fmridetect} }
  \item{anatomic}{ overlay of same dimension as the functional data }
  \item{pos}{ voxel to be marked on output }
  \item{type}{ string. "slice" for slicewise view and "3d" for
    3d view.}
  \item{device}{ output device if \code{type} is \code{slice}. "png",
    "jpeg", "ppm", default is "X11"}
  \item{file}{ name of output file if \code{device} is not "X11" }
}
\details{
  Gives sliceswise view of \code{fmridetect} objects with overlay (if
  needed). 3 dimensional view still experimental.
}
\value{
No value is returned
}
\author{ Karsten Tabelow tabelow@wias-berlin.de }
\note{ 3 dimensional plotting requires \code{misc3d} and \code{tkrplot}
  packages. This is still experimental so do not expect overwhelming results.
}
\seealso{ \code{\link{help}} }
\examples{
\dontrun{plot(detect)}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
