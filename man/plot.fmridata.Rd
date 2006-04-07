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
  \item{x}{ object of class "fmridetect" }
  \item{anatomic}{ overlay of same dimension as the functional data }
  \item{maxpvalue}{ maximum p-value for thresholding }
  \item{spm}{ logical. if class is "fmrispm" decide whether to plot
    the t-statistics of the contrast (\code{spm=TRUE}) or the contrast
    parameter itself (\code{spm=FALSE}).}
  \item{pos}{ voxel to be marked on output }
  \item{type}{ string. "slice" for slicewise view and "3d" for
    3d view.}
  \item{device}{ output device if \code{type} is \code{slice}. "png",
    "jpeg", "ppm", default is "X11"}
  \item{file}{ name of output file if \code{device} is not "X11" }
  \item{...}{ additional arguments for plot}
}
\details{
  Gives sliceswise view of "fmridata" objects with anatomic overlay (if
  appropriate, that is for class "fmripvalue"). For objects of class
  "fmrispm" it plots the t-statistics of the contrast if \code{spm} is
  \code{TRUE}, or the contrast parameter estimate otherwise. For objects of class
  "fmridata" only a plot of the data slices itself is produced. If \code{device} is
  specified as  "png", "jpeg", "ppm" output is done to a file. A scale will
  be plotted in the remaining space of the plot in grey/color scale.

  If \code{type} is "3d" a 3 dimensional interactive view opens. Sliders
  to move in the data cube are given ("x", "y", "z", and "t" if class is
  "fmridata" only). Time series are shown if available. For objects
  of class "fmrispm" a slider is created to remove values below a
  cut-off value from the plot. This removes information for voxels with
  smaller signals than cut-off from the plot. Use pvalues for statistical evaluation. If \code{spm} is
  \code{TRUE} the estimated BOLD response together with a confidence
  interval corresponding to \code{maxpvalue} is drawn. For objects of class
  "fmripvalue" the pvalues with overlay are shown.
}
\value{
  If 'type' is "3d" the Tk-object is returned. (Remove the diplay with tkdestroy(object))
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{ 3 dimensional plotting requires \code{tkrplot}
  package.
}
\seealso{ \code{\link{fmri.pvalue}} }
\examples{
\dontrun{plot(pvalue)}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ iplot }% at least one, from doc/KEYWORDS
