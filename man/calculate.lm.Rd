\name{calculate.lm}
\alias{calculate.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model For FMRI Data }
\description{
  Estimate the parameters and variances of a linear model. 
}
\usage{
  calculate.lm(ttt, z, actype = "smooth", hmax = 3.52, vtype = "var",
                 step = 0.01, contrast = c(1), vvector = c(1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttt}{ 4 dimensional data cube, dimension \eqn{(d_x,d_y,d_z,T)}. }
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
}
\details{
  This function performs parameter estimation in the linear model.
  It implements a two step procedure. In the first step residuals
  are obtained and if \code{actype \%in\% c("ac", "accalc", "smooth")} 
  in each voxel an AR(1) model is fitted to 
  the the time series of residuals, and the AR-coefficient is corrected for bias. 
  If \code{actype=="smooth"}
  the AR-coefficients are spatially smoothed using bandwidth \code{hmax}.
  If \code{actype\%in\%c("ac","smooth")} the linear model is prewithened
  using the estimated (smoothed) AR-coefficients. Parameter 
  and variance estimates are then obtained from the prewithened data.
}
\value{
  \item{beta}{Estimated parameters }
  \item{cbeta}{Estimated contrast of parameters}
  \item{var}{Estimated variance of the contrast of parameters.}
  \item{varm}{Covariance matrix of the parameters given by \code{vvector}}
  \item{res}{Residuals of the estimated linear model}
  \item{arfactor}{Estimated autocorrelation parameter}
}
\references{ Worsley, K.J. (2005). Spatial smoothing of autocorrelations
  to control the degrees of freedom in fMRI analysis. NeuroImage,
  26:635-641.

  Worsley, K.J., Liao, C., Aston, J., Petre, V., Duncan,
  G.H., Morales, F., Evans, A.C. (2002). A general statistical analysis
  for fMRI data. NeuroImage, 15:1:15.
}
\author{ Karsten Tabelow (tabelow@wias-berlin.de) }
\seealso{ \code{\link{create.designmatrix}}, \code{\link{create.stimulus}} }
\examples{
  # Example 1
  ttt <- array(rnorm(32*32*32*107),c(32,32,32,107))
  hrf <- create.stimulus(107, c(18, 48, 78), 15, 2)
  z <- create.designmatrix(hrf,2)
  model <- calculate.lm(ttt,z)
  plot(ttt[16,16,16,])
  lines(ttt[16,16,16,] - model$res[16,16,16,],col=2)
}
\keyword{ regression }% at least one, from doc/KEYWORDS

