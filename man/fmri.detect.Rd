\name{fmri.detect}
\alias{fmri.detect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Thresholding }
\description{
  Performs thresholding.
}
\usage{
fmri.detect(spm, pvalue = 0.05, mode = "plog")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ \code{fmrispm} object }
  \item{pvalue}{ maximal p-value for thresholding (default 0.05) }
  \item{mode}{ string describing kind of signal for detected
    voxels. only "plog" implemented. }
}
\details{
  If only a contrast is given in \code{spm}, we use simply t-statistic
  and thresholding using random field theory for the resulting gaussian
  field (sufficiently large number of df - see ref.). If \code{spm} is a
  vector for each voxel, a chisq field is calculated and evaluated (see
  worsley 2006)
}
\value{
  Object of class \code{fmridetect}
  \item{pvlog}{ log of p-value for activated voxels, zero for non
    activated voxels. use with \code{\link{plot}}.}
}
\references{
  Tabelow, K., Polzehl, J. and Spokoiny, V. (2005). \emph{Analysing {fMRI} experiments with structure adaptive smoothing procedures},
  WIAS-Preprint 1079.

  Worsley, K.J., and Taylor, J.E., \emph{Detecting fMRI activation
    allowing for unknown latency of the hemodynamic response},
  NeuroImage 29:649-654 (2006).
}
\author{ Karsten Tabelow tabelow@wias-berlin.de }
\note{Unxepected side effects may occur if spm does meet the
  requirements, especially if vector for every voxel is used, not only a
  contrast.}
\seealso{ \code{\link{fmri.smooth}}, \code{\link{plot.fmridetect}}  }
\examples{
\dontrun{fmri.detect(smoothresult)}
}
\keyword{ htest }% at least one, from doc/KEYWORDS
