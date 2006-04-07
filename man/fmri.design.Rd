\HeaderA{fmri.lm}{Linear Model For FMRI Data}{fmri.lm}
\keyword{regression}{fmri.lm}
\begin{Description}\relax
Estimate the parameters and variances in a linear model.
\end{Description}
\begin{Usage}
\begin{verbatim}
  fmri.lm(data, z, actype = "accalc", hmax = 3.52, vtype = "var",
          step = 0.01, contrast = c(1), vvector = c(1),
          keep = "essential")
\end{verbatim}
\end{Usage}
\begin{Arguments}
\begin{ldescription}
\item[\code{data}] object of class "fmridata" 
\item[\code{z}] designmatrix specifying the expected BOLD response(s) 
and additional components for trend and other effects. 
\item[\code{actype}] string describing the type of handling autocorrelation
of time series. "nonac", "ac", "accalc", "smooth"
\item[\code{hmax}] bandwidth for smoothing autocorrelation parameter if
actype = "smooth" 
\item[\code{vtype}] method of estimating residual variance (only
"var" implemented) 
\item[\code{step}] step size for binning autocorrelations (see
details) 
\item[\code{contrast}] contrast vector 
\item[\code{vvector}] vector defining the parameters for which the
covariance matrix is returned as well as the corresponding length of
the vector \code{cbeta} in each voxel
\item[\code{keep}] string describing the amount of data returned:
"essential", "diagnostic", "all" 
\end{ldescription}
\end{Arguments}
\begin{Details}\relax
This function performs parameter estimation in the linear model.
It implements a two step procedure. After primary estimation of the
parameters in the first step residuals
are obtained. If \code{actype} \code{\%in\%} \code{c("ac", "accalc", "smooth")} 
an AR(1) model is fitted, in each voxel, to 
the time series of residuals. The estimated AR-coefficient is corrected for bias. 
If \code{actype=="smooth"}
the estimated AR-coefficients are spatially smoothed using bandwidth \code{hmax}.
If \code{actype} \code{\%in\%} \code{c("ac", "smooth")} the linear model is prewithened
using the estimated (smoothed) AR-coefficients. Parameter 
and variance estimates are then obtained from the prewithened
data. The argument \code{keep} describes the amount of data which is
returned. If "essential" only the estimated effects 
\deqn{\tilde{\gamma}_i = C^T\tilde{\beta}_i}{}
and their
estimated variances are returned. "all" gives the full data, including
residuals, temporal autocorrelation.
If \code{vvector} is given and has length greater than 1, the
covariance matrix for the stimuli given therein are returned
(\code{varm}) and \code{vwghts} contains an estimate for the ratio of
the variances of the parameter for the stimuli indicated in
\code{vvector}. \code{cbeta} then contains the corresponding parameter
estimates and thus is a vector of corresponding length in each voxel.
\end{Details}
\begin{Value}
object with class attributes "fmrispm" and "fmridata"
\begin{ldescription}
\item[\code{beta}] Estimated parameters 
\item[\code{cbeta}] Estimated contrast of parameters
\item[\code{var}] Estimated variance of the contrast of parameters.
\item[\code{varm}] Covariance matrix of the parameters given by \code{vvector}
\item[\code{res}] Residuals of the estimated linear model
\item[\code{arfactor}] Estimated autocorrelation parameter
\item[\code{scorr}] spatial correlation of data 
\item[\code{weights}] ratio of voxel dimensions 
\item[\code{vwghts}] ratio of estimated variances for the stimululi given by
\code{vvector}
\item[\code{rxyz}] array of smoothness from estimated correlation for each
voxel in resel space (for analysis without smoothing)
\item[\code{hrf}] Expected BOLD response for contrast 
\end{ldescription}
\end{Value}
\begin{Note}\relax
\code{vvector} is intended to be used for delay of the HRF using
its first derivative. Do not mix with the \code{contrast} argument, since
unexpected side effects may occur. Look out for updates of this package.
\end{Note}
\begin{Author}\relax
Karsten Tabelow \email{tabelow@wias-berlin.de}
\end{Author}
\begin{References}\relax
Worsley, K.J. (2005). Spatial smoothing of autocorrelations
to control the degrees of freedom in fMRI analysis. NeuroImage,
26:635-641.

Worsley, K.J., Liao, C., Aston, J., Petre, V., Duncan,
G.H., Morales, F., Evans, A.C. (2002). A general statistical analysis
for fMRI data. NeuroImage, 15:1-15.
\end{References}
\begin{SeeAlso}\relax
\code{\LinkA{fmri.design}{fmri.design}}, \code{\LinkA{fmri.stimulus}{fmri.stimulus}}
\end{SeeAlso}
\begin{Examples}
\begin{ExampleCode}
  # Example 1
  data <- list(ttt=array(rnorm(32*32*32*107),c(32,32,32,107)),
               mask=array(1,c(32,32,32)))
  class(data) <- "fmri.data"
  hrf <- fmri.stimulus(107, c(18, 48, 78), 15, 2)
  z <- fmri.design(hrf,2)
  model <- fmri.lm(data,z,keep="all")
  plot(data$ttt[16,16,16,])
  lines(data$ttt[16,16,16,] - model$res[16,16,16,],col=2)
\end{ExampleCode}
\end{Examples}

