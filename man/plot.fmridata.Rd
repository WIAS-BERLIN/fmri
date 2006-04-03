\name{plot.fmridata}
\alias{plot.fmridata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  Plot the fMRI detection result and thresholding.
}
\usage{
plot.fmridata(x, anatomic = NULL, maxpvalue = 0.05,
              spm = TRUE, pos = c(-1, -1, -1), type = "slice",
              device = "X11", file = "plot.png",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{fmridetect} }
  \item{anatomic}{ overlay of same dimension as the functional data }
  \item{maxpvalue}{ maximum p-value for thresholding }
  \item{spm}{ logical. if class is '"fmrispm"' decide whether to plot
    the t-statistics of the contrast or the contrast parameter itself.}
  \item{pos}{ voxel to be marked on output }
  \item{type}{ string. "slice" for slicewise view and "3d" for
    3d view.}
  \item{device}{ output device if \code{type} is \code{slice}. "png",
    "jpeg", "ppm", default is "X11"}
  \item{file}{ name of output file if \code{device} is not "X11" }
  \item{...}{ additional arguments for plot}
}
\details{
  Gives sliceswise view of '"fmridata"' objects with anatomic overlay (if
  appropriate, that is for class '"fmripvalue"'). Objects of class
  '"fmrispm"' plot the t-statistics of the contrast if \code{spm} is
  \code{TRUE}, or the contrast parameter estimate otherwise. Objects of class
  '"fmridata"' only plot the data slices itself. If \code{device} is
  given as  "png", "jpeg", "ppm" output is done to a file. A scale will
  be plotted in the remaining space of the plot.

  If \code{type} is "3d" a 3 dimensional interactive view opens. Sliders
  to move in the data cube are given ("x", "y", "z", and "t" if class is
  '"fmridata"' only). Time series are shown if available. For objects
  of class '"fmrispm"' a threshold slider is created to remove below
  threshold values from the plot. This is to look at data properties
  only! Use pvalues for statistical evaluation. If \code{spm} is
  \code{TRUE} the estimated BOLD response together with a confidence
  interval corresponding to \code{maxpvalue} is drawn. For Objects of class
  '"fmripvalue"' the pvalues with overlay are shown.
}
\value{
  If 'type' is '"3d"' the Tk-object is returned. (Destroy with tkdestroy(object))
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{ 3 dimensional plotting requires \code{tkrplot}
  packages.
}
\seealso{ \code{\link{fmri.pvalue}} }
\examples{
\dontrun{plot(pvalue)}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ iplot }% at least one, from doc/KEYWORDS
