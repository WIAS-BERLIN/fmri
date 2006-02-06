\name{fmriplot.signal}
\alias{fmriplot.signal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ FMRI Results }
\description{
  Plot fMRI results with underlay slicewise.
}
\usage{
fmriplot.signal(signal, anatomic = array(0, dim = dim(signal)), pos = c(-1, -1, -1), device = "X11", file = "plot.png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal}{ 3 dimensional array of signal strengths (or p-value or
  ...). }
  \item{anatomic}{ Underlay for signal. (anatomic) }
  \item{pos}{ Mark this voxel }
  \item{device}{ Output device. "X11", "png", "jpeg", "ppm" }
  \item{file}{ Name of the file for output if device is not X11 }
}
\details{
  Plots slicewise the results given by \code{signal} with underlay given
  in \code{anatomic}. Voxel \code{pos} is marked.

  Output can be X11 and a graphics file \code{device}.
}
\value{
  Nothing is returned}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{fmriplot.timeseries}}, \code{\link{fmriplot.slices3d}}, \code{\link{fmriplot.diff}} }
\examples{
  \dontrun{fmriplot.signal(signal, anatomic = array(0, dim = dim(signal)), pos = c(-1, -1, -1), device = "X11", file = "plot.png")}
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
