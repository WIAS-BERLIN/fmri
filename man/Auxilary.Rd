\name{Auxilary}
\alias{get3Dh.gauss}
\alias{Spatialvar.gauss}
\alias{get.bw.gauss}
\alias{Varcor.gauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Auxilary functions} 
\description{
  Auxilary functions for caculation of variance reduction, correlation,
  bandwiths corresponding to a specified variance reduction or to a specified correlation.
}
\usage{
get3Dh.gauss(vred, h0, vwghts, step = 1.002, interv = 1)
}
\usage{
Spatialvar.gauss(h, h0, d, interv = 1)
}
\usage{
get.bw.gauss(corr, step = 1.001, interv = 2)
}
\usage{
Varcor.gauss(h, interv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vred}{ variance reduction }
  \item{h}{ vector (length 1 or d) of bandwidths in FWHM units }
  \item{h0}{ vector (length 1 or d) of bandwidths in FWHM units }
  \item{vwghts}{ vector of weights }
  \item{step}{ step size for iteration }
  \item{d}{ dimension of the data (1, 2 or 3) }
  \item{interv}{ number of intervals per unit for discretization of the Gaussian kernel. 
    \code{interv=1} corresponds to kernel smoothing, \code{interv>1} allows to handle
    situations where Gaussian smoothing preceeds discretisation into voxels. The difference is
    relevant for small bandwidths. }
}
\details{ \item{get3Dh.gauss}{ Calculate the bandwidth corresponding to the amount of variance
  reduction achieved locally.}
        \item{get.bw.gauss} {Calculate the   bandwidth for the Gaussian kernel corresponding to a given correlation.
	   Use \code{interv=1} if the correlation results from kernel smoothing and \code{interv>1} for
	   intrinsic correlation (smoothing preceeding discretisation into voxel). }
 \item{Varcor.gauss} {Calculates a correction for the variance estimate obtained by
  \eqn{(IQRdiff(y)/1.908)^2}{(IQRdiff(y)/1.908)^2} in case of colored noise that was produced by
  smoothing with gaussian kernel and bandwidth \code{h}.}
}
\value{
  bandwidth
}
}
\references{  }
\author{Joerg Polzehl, \email{polzehl@wias-berlin.de}}
\note{These function sare for internal use only.
}
\seealso{ \code{\link{vaws3D}} }
\examples{
# get bandwidths (h_x,h_y,h_z) of a Gaussian kernel applied to spatially correlated data that achieves a variance reduction 
# (with respect to  independent data) by factor vred  
# 
\dontrun{get3Dh.gauss(vred=c(.04,.08,.1), h0=c(.5,.7,.3), vwghts=c(1,1,1), step = 1.002, interv = 1)}
# get variance reduction (with respect to  independent data) achieved by a Gaussian kernel with bandwidths h=(h_x,h_y,h_z)  
# applied to spatially correlated data 
# 
\dontrun{Spatialvar.gauss(h=c(2.5,2.5,1.5), h0=c(.5,.7,.3), d=3, interv = 1)}
# get bandwidth of a Gaussian kernel (1D) that produces thespecified spatial correlation
# 
\dontrun{get.bw.gauss(corr=.4, step = 1.001, interv = 2)}
#  calculate the correction factor needed for a variance estimate by
#  (IQRdiff(y)/1.908)^2}{(IQRdiff(y)/1.908)^2}
#  in case of spatial correlation of .4,.5 and .3 in the three coordinate directions
\dontrun{
Varcor.gauss(h=c(get.bw.gauss(.4),get.bw.gauss(.5),get.bw.gauss(.3)), interv = 1)
}
}
\keyword{internal}% at least one, from doc/KEYWORDS
