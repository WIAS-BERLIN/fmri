\name{fmriplot.timeseries}
\alias{fmriplot.timeseries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ FMRI Results }
\description{
  Plot fMRI time series and fitted model.
} 
\usage{
fmriplot.timeseries(fmri, fit, sigvec, tvec, thr, origin, 
                       device = "X11", file = "plot.png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fmri}{ 2+1 dimensional fMRI data }
  \item{fit}{ 2+1 dimensional data of the fitted model }
  \item{sigvec}{ 2 dimensional array with values of the signal }
  \item{tvec}{ Value of the t-statistic }
  \item{thr}{ Threshold value }
  \item{origin}{ Lower left voxel coordinate of the \code{fmri} data }
  \item{device}{ Output device. "X11", "png", "jpeg", "ppm" }
  \item{file}{ Name of the file for output if device is not X11 }
}
\details{
  This function shows a 2 dimensional grid of time series together with
  the fitted lines of the model for the time series. Voxel coordinates
  and estimated signal strength together with value of the t-statistics
  are plotted. Signifikant voxels are indicated by red dots. Threshold
  to be used (globally) is given by \code{threshold}
}
\value{
  Nothing is returned.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{ This function is very spezialized to view the output of fmri
  analysis. It may not be of much help.}
\seealso{ \code{\link{fmriplot.signal}}, \code{\link{fmriplot.slices3d}}, \code{\link{fmriplot.diff}}  }
\examples{
  \dontrun{fmriplot.timeseries(fmri, fit, sigvec, tvec, thr, origin, 
                   device = "X11", file = "plot.png")}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
