\name{resel}
\alias{resel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Expected Euler Characteristic densities }
\description{
  Calculate voxel size in resel space.
}
\usage{
resel(voxeld, hmax, hv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{voxeld}{ voxeldimension in mm }
  \item{hmax}{ FWHM bandwidth of used gaussian kernel for smoothing }
  \item{hv}{ correction factor (not needed) }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  Calculate voxel size in resel space. See references for detailed
  description. 
}
\references{ K.J. Worsley, S. Marrett, P. Neelin, A.C. Vandal,
  K.J. Friston, and A.C. Evans (1996), A Unified Statistical Approach
  for Determining Significant Signals in Images of Cerebral Activation, 
  Human Brain Mapping, 4:58-73.}
\author{ Karsten Tabelow (tabelow@wias-berlin.de }
\seealso{ \code{\link{pvalue}}, \code{\link{threshold}} }
\examples{
  # Example 1
  resel(4,10)
}
\keyword{ misc }% at least one, from doc/KEYWORDS
