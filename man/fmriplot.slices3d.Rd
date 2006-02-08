\name{fmriplot.slices3d}
\alias{fmriplot.slices3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ FMRI Results }
\description{
  Plot 3D slices of the FMRI results.
}
\usage{
fmriplot.slices3d(signal, anatomic = array(0, dim = dim(signal)), 
                                        pos = c(-1, -1, -1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal}{ 3 dimensional array of signal strengths (or p-value or
  ...). }
  \item{anatomic}{ Underlay for signal. (anatomic) }
  \item{pos}{ Mark this voxel }
}
\details{
  Plot 3D slices of the FMRI results with underlay. If \code{pos} is
  negativ, no voxel is marked.
}
\value{
  Nothing is returned}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{
  This code is still experimental! Should rely on tkrplot only in future
  versions.
}
\seealso{ \code{\link{fmriplot.signal}},
  \code{\link{fmriplot.timeseries}}, \code{\link{fmriplot.diff}} }
\examples{
  \dontrun{fmriplot.slices3d(signal, anatomic = array(0, dim = dim(signal)), 
                                  pos = c(-1, -1, -1))}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ dynamic }% __ONLY ONE__ keyword per line
