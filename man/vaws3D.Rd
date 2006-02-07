\name{vaws3D}
\alias{vaws3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Three dimensional adaptive smoothing of vector valued data}
\description{Performs three dimensional adaptive smoothing of vector valued data
using the Propagation-Separation approach from Polzehl \& Spokoiny (2006). 
}
\usage{
vaws3D(y, qlambda = NULL, qtau = NULL, lkern = "Triangle", skern="Exp", aggkern = "Uniform", sigma2 = NULL, hinit = NULL, hincr = NULL, hmax = NULL, lseq = NULL, heta = NULL, u = NULL, graph = FALSE, demo = FALSE, wghts = NULL, spmin = 0, spmax = 5, scorr = 0, vwghts = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Data array of dimension \code{c(dx,dy,dz,dv)}. The
  first three values specify the size of the 3D data volume. \code{dv}
  specifies the length of the data vector in each voxel.  }
  \item{qlambda}{ \code{qlambda} is the main smoothing parameter. It determines the scaling parameter \code{\lambda} in the statistical penalty
  (see description of the algorithm in the references) as \code{qchisq(qlambda,sum(vwghts^2)^2/sum(vwghts^4))}. \code{qlambda=1} disables adaptation, i.e. the 
  resulting estimate is a kernel estimate using the largest inspected bandwidth less or equal to \code{hmax}.
  The default value for \code{qlambda} is selected as the smallest value that fulfils a propagation condition.
  This means that in a parametric, in this case constant, model the resulting estimate is, for a large \code{hmax} almost parametric. 
   Larger values of 
  \code{qlambda} lead to less sensitivity to structural differences, smaller values may lead to a random structure (segmentation) of the resulting estimate.}
  \item{qtau}{Stagewise Aggregation, see Belomestny and Spokoiny (2004) is used as an adaptive control step if \code{qtau<1}. 
   \code{qtau} determines the scaling parameter \code{\tau} in the Stagewise Aggregation algorithm in the same way as \code{qlambda} does
   for the PS-approach. Default values are again selected by a propagation condition.}
  \item{lkern}{ \code{lkern} specifies the location kernel. Defaults to "Triangle", other choices are "Gaussian", "Quadratic", "Cubic" and "Uniform".
    Note that the location kernel is applied to \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression.}
  \item{skern}{ \code{lkern} specifies the kernel for the statistikal penalty. Defaults to "Exp", the alternative is "Triangle".
    \code{lkern="Triangle"} allows for much faster computation (saves up to 50%) but needs some additional testing.}
  \item{aggkern}{\code{aggkern} specifies the kernel for the statistical panalty in stagewise aggregation. Defaults to "Uniform", the alternative choice is
   "Triangle" }
  \item{sigma2}{\code{sigma2} allows to specify the variance \eqn{\sigma^2_i} of the first component of the data vector in each voxel \eqn{i}. 
   Defaults to \code{NULL}. In this case a homoskedastic variance estimate is generated. If \code{length(sigma2)==1} a homoskedastic
   variance is assumed. }
  \item{hinit}{\code{hinit} specifies the initial bandwidth. Defaults to \code{hinit=1} (\code{hinit=1/4} in case of \code{lkern="Gaussian"}).  }
  \item{hincr}{\code{hincr} specifies the factor used to increase the size of local neigborhoods after each iteration. The bandwidth is increased by
   a factor \code{hincr^(1/3)}. Defaults to \code{hincr=1.25}   }
  \item{hmax}{\code{hmax} specifies the maximal bandwidth. Defaults to \code{hmax= 5}}
  \item{lseq}{\code{lseq} allows to increase the value of the scaling parameter \code{\lambda} for the first \code{length(lseq)} iteration steps by
  the factor specified in \code{lseq}. Defaults to NULL. In this case a default sequence is used that fulfils the propagation condition.}
  \item{heta}{\code{heta} specifies the minimal bandwidth to use with stagewise aggregation. }
  \item{u}{\code{u} allows to specify the true parameter. This is only used to test the algorithm and to select the smoothing parameters 
    \code{qlambda} and \code{qtau} by a propagation condition. If \code{u} is specified MSE and MAE of the estimates are 
    printed for each iteration step. }
  \item{graph}{If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}.  }
  \item{demo}{If \code{demo=TRUE} the function pauses after each iteration. Defaults to \code{demo=FALSE}. }
  \item{wghts}{\code{wghts} specifies the  diagonal elements of a weight matrix to adjust for different distances between grid-points
  in different coordinate directions, i.e. allows to define a more appropriate metric in the design space. }
  \item{spmin}{\code{spmin} specifies parameter value in the kernel \eqn{K_s(x)=\min(0,e^{-spmax/(spmax-spmin)(x-spmin)}I_{x<spmax}}. Defaults
  to \code{spmin=0}.}
  \item{spmax}{\code{spmax} specifies the cut point in the kernel \eqn{K_s(x)=\min(0,e^{-spmax/(spmax-spmin)(x-spmin)}I_{x<spmax}} }. Defaults
  to \code{spmax=5}.}
  \item{scorr}{\code{scorr} specifies the spatial correlation in the three coordinate directions. Defaults to \code{scorr=c(0,0,0)}. }
  \item{vwghts}{\code{vwghts} specifies weights for the components of the data vector. These weights are used in the evaluation of the statistical
  penalty \eqn{s_ij=N_i/\lambda/\sigma^2_i sum_k vwghts[k] (theta_{ik}-\theta_{jk})^2}{} where \eqn{k} denotes the \eqn{k}-th component of the 
  vector. The parameter can be used to adjust for different variability of the components of the data vector. Defaults to \code{vwghts=rep(1,dv)}.  }
}
\details{The function implements the propagation separation approach to 
nonparametric smoothing (formerly introduced as Adaptive weights smoothing) 
for multivariate regression models regression with additive "Gaussian" errors on a 3D grid. A homoskedastic 
or heteroskedastic model is used depending on the content of \code{sigma2}. 
\code{qlambda>=1} provides the stagewise aggregation procedure from Belomestny and Spokoiny (2004).
\code{qtau>=1} provides Adaptive weights smoothing without control by stagewise aggregation.  

The essential parameter in the procedure is \code{qlambda}. This parameter has an 
   interpretation as a significance level of a test for equivalence of two local
   parameter estimates.
   Default values provided are choosen to fulfil the propagation, i.e. in case of a 
   constant (global) parameter value and large \code{hmax} the procedure should 
   provide, with a high probability, the global (parametric) estimate.
   More formally we require the parameter \code{qlambda} and eventually \code{lseq}
   to be specified such that
   \code{\bf{E} |\hat{\theta}^k - \theta| \le (1+\alpha) \bf{E} |\tilde{\theta}^k - \theta|}
   where \code{\hat{\theta}^k} is the aws-estimate in step \code{k} and \code{\tilde{\theta}^k}
   is corresponding nonadaptive estimate using the same bandwidth (\code{qlambda=1}).
   Default values are selected to fulfil this condition for \code{\alpha=0.1}.
   
   The optimal values only slightly depend on the model parameters, i.e. the
   default parameters should work in most situations. Larger values of \code{qlambda}
   may lead to oversmoothing, small values of \code{qlambda} lead to a random segmentation
   of homogeneous regions. 
   
   The numerical complexity of the procedure is mainly determined by \code{hmax}. The number
   of iterations is \code{d*log(hmax/hinit)/log(hincr)} with \code{d} being the dimension 
   of \code{y}. Comlexity in each iteration step is \code{Const*hakt*n} with \code{hakt}
   being the actual bandwith in the iteration step and \code{n} the number of design points.
   \code{hmax} determines the maximal possible variance reduction.
}
\value{
  \item{theta}{Array of dimension \code{dim(y)} containing the smoothed data. }
  \item{ni}{Array of dimension \code{dim(y)[1:3]} containing the voxelwise sum of
  weights.}
  \item{qi}{Array of dimension \code{dim(y)[1:3]} containing the voxelwise sum of
  squared weights.}
  \item{cgh}{The value \code{cgh}.}
  \item{var}{The estimated variances of \code{theta[,,,1]}.}
  \item{vred}{The amount of variance reduction achieved in each voxel.}
  \item{y}{The data.}
  \item{hmax}{The largest bandwidth used.}
  \item{mae}{Mean absolute error if \code{!is.null(u)}.}
  \item{mse}{Mean squared error if \code{!is.null(u)}.}
  \item{lseq}{The sequence specified in \code{lseq}.}
  \item{args}{The arguments used in the call of function vaws3D}  
}
\references{ 
\item{ }{Polzehl, J. and Spokoiny, V. (2004a). \emph{Propagation-Separation Approach for Local Likelihood Estimation}, 
Probab. Theory & Relat. Fields, in print.}
\item{ }{Belomestny, D. and Spokoiny, V. (2004a). \emph{Local likelihood modeling via stagewise aggregation}, 
WIAS-Preprint 1000.}
\item{ }{Tabelow, K., Polzehl, J. and Spokoiny, V. (2005). \emph{Analysing {fMRI} experiments with structure adaptive smoothing procedures},
WIAS-Preprint 1079.} 
 }
\author{Joerg Polzehl, \email{polzehl@wias-berlin.de}, 
\url{http://www.wias-berlin.de/project-areas/stat/projects/adaptive-image-processing.html} }
#\note{} 
}
\seealso{\code{\link{perform.aws}}}
\examples{
#
#  create an 3D phantom
#
mask<-matrix(0,65,65)
mask[5:10,5:10]<-1
mask[7:8,7:8]<-0
mask[8:10,8:10]<-0
mask[14:17,14:17]<-1
mask[16:17,16:17]<-0
mask[21:23,21:23]<-1
mask[22:23,23]<-0
mask[23,22]<-0
mask[27:28,27:28]<-1
mask[28,28]<-0
mask[5:7,29:33]<-1
mask[7,32:33]<-0
mask[14:15,30:33]<-1
mask[15,30]<-0
mask[21,31:33]<-1
mask[22,33]<-1
mask[27,32:33]<-1
mask[29:33,5:7]<-1
mask[32:33,7]<-0
mask[30:33,14:15]<-1
mask[30,15]<-0
mask[31:33,21]<-1
mask[33,22]<-1
mask[32:33,27]<-1
mask[34:65,1:33]<-mask[32:1,1:33]
mask[1:33,34:65]<-mask[1:33,32:1]
mask[34:65,34:65]<-mask[32:1,32:1]
image(mask)
size<-mask
signal<-1
factor<-1.1
size[29:37,38:65]<-signal
size[38:65,38:65]<-signal*factor
size[38:65,29:37]<-signal*factor^2
size[38:65,1:28]<-signal*factor^3
size[29:37,1:28]<-signal*factor^4
size[1:28,1:28]<-signal*factor^5
size[1:28,29:37]<-signal*factor^6
size[1:28,38:65]<-signal*factor^7
size<-size*mask
data<-array(0,c(65,65,10,1))
for(i in 1:10) data[,,i,1]<-size
mdata<-2*data
#  
#   add some noise
#
data<-mdata+array(rnorm(data),dim(data))
#
#   smooth the data ( takes about 6s on a Pentium IV 3.2GHz )
#
ttt<-vaws3D(data,sigma2=array(1,dim(data)[-4]),skern="Triangle",hmax=4,graph=FALSE,u=mdata,qtau=1)$theta
}
\keyword{smooth}% at least one, from doc/KEYWORDS
