\name{fmri.lm}
\alias{fmri.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model For FMRI Data }
\description{
  Estimate the parameters and variances of a linear model. 
}
\usage{
  fmri.lm(data, z, actype = "accalc", hmax = 3.52, vtype = "var",
          step = 0.01, contrast = c(1), vvector = c(1),
          keep = "essential")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ object of class \code{fmri.data} }
  \item{z}{ designmatrix specifying the expected hemodynamic response 
  and aadditional components for trend and other effects. Dimension \eqn{(p,T)}. }
  \item{actype}{ string describing the type of handling autocorrelation
    of time series. "nonac", "ac", "accalc", "smooth"}
  \item{hmax}{ bandwidth for smoothing autocorrelation parameter if
    actype = "smooth" }
  \item{vtype}{ type of estimation of estimating residual variance (only
  "var" implemented) }
  \item{step}{ step size for binning autocorrelation handling (see
    details) }
  \item{contrast}{ contrast vector }
  \item{vvector}{ vector defining the parameters for which the
    covariance matrix is returned. }
  \item{keep}{ string describing the amount of data returned:
    "essential", "diagnostic", "all" }
}
\details{
  This function performs parameter estimation in the linear model.
  It implements a two step procedure. In the first step residuals
  are obtained and if \code{actype} \code{\%in\%} \code{c("ac", "accalc", "smooth")} 
  in each voxel an AR(1) model is fitted to 
  the the time series of residuals, and the AR-coefficient is corrected for bias. 
  If \code{actype=="smooth"}
  the AR-coefficients are spatially smoothed using bandwidth \code{hmax}.
  If \code{actype\%in\%c("ac","smooth")} the linear model is prewithened
  using the estimated (smoothed) AR-coefficients. Parameter 
  and variance estimates are then obtained from the prewithened
  data. The argument \code{keep} describes the amount of data which is
  returned. If "essential" only the contrast of parameters and its
  estimated variance are returned. "all" gives the full data, including
  residuals, temporal autocorrelation.
  If \code{vvector} is given and has length greater than 1, the
  covariance matrix for the stimuli given therein are returned
  (\code{varm}) and \code{vwghts} contains an estimate for the ratio of
  the variances of the parameter for the stimuli indicated in \code{vvector}.
}
\value{
  object of class \code{fmrispm}
  \item{beta}{Estimated parameters }
  \item{cbeta}{Estimated contrast of parameters}
  \item{var}{Estimated variance of the contrast of parameters.}
  \item{varm}{Covariance matrix of the parameters given by \code{vvector}}
  \item{res}{Residuals of the estimated linear model}
  \item{arfactor}{Estimated autocorrelation parameter}
  \item{scorr}{ spatial correlation of data }
  \item{weights}{ ratio of voxel dimensions }
  \item{vwghts}{ ratio of estimated variances for the stimululi given by
  \code{vvector}}
  \item{rxyz}{array of resel counts from estimated correlation for each
    voxel (for analysis without smoothing)}
  \item{smooth}{ FALSE }
}
\references{ Worsley, K.J. (2005). Spatial smoothing of autocorrelations
  to control the degrees of freedom in fMRI analysis. NeuroImage,
  26:635-641.

  Worsley, K.J., Liao, C., Aston, J., Petre, V., Duncan,
  G.H., Morales, F., Evans, A.C. (2002). A general statistical analysis
  for fMRI data. NeuroImage, 15:1:15.
}
\author{ Karsten Tabelow (tabelow@wias-berlin.de) }
\note{ \code{vvector} is intended to be used for delay of HRF using
  first derivative. Do not mix with \code{contrast} argument, since
  unexpected side effects may occur. Look out for updates of this package.
}
\seealso{ \code{\link{fmri.design}}, \code{\link{fmri.stimulus}} }
\examples{
  # Example 1
  data <- list(ttt=array(rnorm(32*32*32*107),c(32,32,32,107)),mask=array(1,c(32,32,32)))
  class(data) <- "fmri.data"
  hrf <- fmri.stimulus(107, c(18, 48, 78), 15, 2)
  z <- fmri.design(hrf,2)
  model <- fmri.lm(data,z,keep="all")
  plot(data$ttt[16,16,16,])
  lines(data$ttt[16,16,16,] - model$res[16,16,16,],col=2)
}
\keyword{ regression }% at least one, from doc/KEYWORDS
