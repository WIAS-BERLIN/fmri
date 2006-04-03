\name{fmri.pvalue}
\alias{fmri.pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ P-values }
\description{
  Determines p-values.
}
\usage{
fmri.pvalue(spm, mode="basic", delta=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ \code{fmrispm} object }
  \item{mode}{ type of pvalue definition }
  \item{delta}{ physically meaningful range of latency for HRF }
}
\details{
  If only a contrast is given in \code{spm}, we use simply t-statistic
  and thresholding using random field theory for the resulting gaussian
  field (sufficiently large number of df - see ref.). If \code{spm} is a
  vector for each voxel, a chisq field is calculated and evaluated (see
  worsley 2006). If \code{delta} is given, a cone statistics is
  used. The parameter \code{mode} allows for different kinds of p-value
  calculation. "basic" corresponds to a global definition of the
  resel counts based on the amount of smoothness achieved by a equivalent
  Gaussian filter. The propagation condition ensures, that under the
  hypothesis \deqn{\hat{\Theta} = 0} the adaptive smoothing perform like a
  non adaptive filter with the same kernel function. "local"
  corresponds to a more conservative approach, where the p-value is
  derived from the estimated local smoothness that has been achieved by
  the adaptive smoothing. "global" takes a global mean.  }
\value{
  Object of class \code{fmripvalue} and \code{fmridata}
  \item{pvalue}{ p-value. use with \code{\link{plot}} for thresholding.}
  \item{weights}{ voxelsize ratio }
  \item{dim}{ data dimension }
  \item{hrf}{ expected BOLD response for contrast (single stimulus only)}
}
\references{
  Tabelow, K., Polzehl, J. and Spokoiny, V. (2005). \emph{Analysing {fMRI}
    experiments with structure adaptive smoothing procedures},
  WIAS-Preprint 1079. 
  
  Worsley, K.J., and Taylor, J.E., \emph{Detecting fMRI activation
    allowing for unknown latency of the hemodynamic response},
  NeuroImage 29:649-654 (2006).
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{Unxepected side effects may occur if spm does meet the
  requirements, especially if vector for every voxel is used, not only a
  contrast.}
\seealso{ \code{\link{fmri.smooth}}, \code{\link{plot.fmridata}}  }
\examples{
\dontrun{fmri.pvalue(smoothresult)}
}
\keyword{ htest }% at least one, from doc/KEYWORDS
