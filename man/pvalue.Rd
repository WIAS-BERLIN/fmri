\name{pvalue}
\alias{pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Expected Euler Characteristic densities }
\description{
  Calculate the Expected Euler Characteristic densities for a 3d gaussian,
  t-, or chisq-field.
}
\usage{
pvalue(z, i, j, k, rx, ry, rz, type = "norm", df = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ value of the field }
  \item{i}{ x-dimension of the field }
  \item{j}{ y-dimension of the field }
  \item{k}{ z-dimension of the field }
  \item{rx}{ voxel size in resel space x-direction }
  \item{ry}{ voxel size in resel space y-direction }
  \item{rz}{ voxel size in resel space z-direction }
  \item{type}{ string for the type of field, currently implemented:
    "norm", "t", "chisq" }
  \item{df}{ number of degrees of freedom for t- or chisq-field }
}
\details{
  The function calculates the Expected Euler Characteristic densities
  for a 3d gaussian, t-, or chisq-field. If lower dimensional field is
  needed simply set k=1 and j=1.

  Boxes are considered as resel shapes.

  See references for detailed description and formulas. 
}
\value{
  p-value for the field value z (if sufficiently large) according to
  Expected Euler Characteristic density.
}
\references{ K.J. Worsley, S. Marrett, P. Neelin, A.C. Vandal,
  K.J. Friston, and A.C. Evans (1996), A Unified Statistical Approach
  for Determining Significant Signals in Images of Cerebral Activation, 
  Human Brain Mapping, 4:58-73.}
\author{ Karsten Tabelow (tabelow@wias-berlin.de }
\seealso{ \code{\link{threshold}},\code{\link{resel}}  }
\examples{
  # Example 1
  z <- seq(0,10,0.1)
  rx <- resel(4,10) # voxel dimension 4mmm in x-direction
  ry <- resel(4,10) # voxel dimension 4mmm in y-direction
  rz <- resel(5,10) # voxel dimension 5mmm in z-direction
  p <- pvalue(z,64,64,30,rx,ry,rz)
  plot(z,p,type="l")
}
\keyword{ distribution }% __ONLY ONE__ keyword per line
