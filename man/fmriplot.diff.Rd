\name{fmriplot.diff}
\alias{fmriplot.diff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ FMRI Results }
\description{
  Plot difference between fMRI results with underlay slicewise.
}
\usage{
fmriplot.diff(mask1, mask2, anatomic = array(0, dim = dim(mask1)), pos = c(-1, -1, -1), device = "X11", file = "plot.png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mask1}{ 3 dimensional array mask of significant results by
    method 1 }
  \item{mask2}{ 3 dimensional array mask of significant results by
    method 2 }
  \item{anatomic}{ Underlay for signal. (anatomic) }
  \item{pos}{ Mark this voxel }
  \item{device}{ Output device. "X11", "png", "jpeg", "ppm" }
  \item{file}{ Name of the file for output if device is not X11 }
}
\details{
  Plots the difference of the results of two fmri activation detection
  methods with underlay slicewise. Voxels found significant by method 1
  only are shown with green colour. Voxels found significant by method 2
  only are shown with blue colour. Voxels found significant by both
  methods are shown with red colour.

  \code{mask1} and \code{mask2} should be given as 3 dimensional logical
  arrays.
} 
\value{
  Nothing is returned
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{fmriplot.signal}}, \code{\link{fmriplot.slices3d}}, \code{\link{fmriplot.timeseries}} }
\examples{
\dontrun{fmriplot.diff(mask1, mask2, anatomic = array(0, dim = dim(mask1)), pos = c(-1, -1, -1), device = "X11", file = "plot.png")}
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
