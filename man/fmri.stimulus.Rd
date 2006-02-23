\name{fmri.stimulus}
\alias{fmri.stimulus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model for FMRI Data }
\description{
  Creates expected BOLD response for a given task indicator function.
}
\usage{
  fmri.stimulus(scans = 1, onsets = c(1), length = 1, rt = 3, mean =
                 TRUE, a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scans}{ number of scans }
  \item{onsets}{ vector of onset times (in scans) }
  \item{length}{ length of ON stimulus (in scans) }
  \item{rt}{ time between scans in seconds (TR) }
  \item{mean}{ logical. if TRUE apply meann correction to resulting vector }
  \item{a1}{ parameter of hemodynamic response function (see details) }
  \item{a2}{ parameter of hemodynamic response function (see details) }
  \item{b1}{ parameter of hemodynamic response function (see details) } 
  \item{b2}{ parameter of hemodynamic response function (see details) }
  \item{cc}{ parameter of hemodynamic response function (see details) }
}
\details{
  The functions calculates and expected BOLD response for the task
  indicator function given as argument as convolution with the
  hemodynamic response function. The latter is modelled difference
  between two gamma functions as given in the reference (with defaults
  given there):
  
  c1 * exp(-(x - d1)/b1) - c2 * exp(-(x - d2)/b2)

  The parameters of this function can be changed in the call of
  fmri.stimulus.

  Dimension of the function value is set to c(scans,1).

  If \code{mean} is TRUE (default) the resulting vector is
  mean-corrected.
}
\value{
  Vector with \code{scans} elements. Dimension is set to c(scans,1).
}
\references{ Worsley, K.J., Liao, C., Aston, J., Petre, V., Duncan,
  G.H., Morales, F., Evans, A.C. (2002). A general statistical analysis
  for fMRI data. NeuroImage, 15:1:15. }
\author{ Karsten Tabelow (tabelow@wias-berlin.de) }
\seealso{ \code{\link{fmri.design}}, \code{\link{fmri.lm}} }
\examples{
  # Example 1
  hrf <- fmri.stimulus(107, c(18, 48, 78), 15, 2)
  z <- fmri.design(hrf,2)
  dim(hrf)

}
\keyword{ regression }% at least one, from doc/KEYWORDS
\keyword{ design }% __ONLY ONE__ keyword per line
