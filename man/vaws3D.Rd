\name{vaws3D}
\alias{vaws3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Three dimensional adaptive smoothing of vector valued data}
\description{Performs three dimensional adaptive smoothing of vector valued data
using the Propagation-Separation approach from Polzehl \& Spokoiny (2006). 
}
\usage{
vaws3D(y, qlambda = NULL, qtau = NULL, lkern = "Triangle", aggkern = "Uniform", sigma2 = NULL, hinit = NULL, hincr = NULL, hmax = NULL, lseq = NULL, heta = NULL, u = NULL, graph = FALSE, demo = FALSE, wghts = NULL, spmin = 0, spmax = 5, scorr = 0, vwghts = 1)
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
  penalty \eqn{s_ij=N_i/\lambda/\sigma^2_i sum_k vwghts[k] (\theta_{ik}-\theta_{jk})^2} where \eqn{k} denotes the \eqn{k}-th component of the 
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
    z <- list(theta = tobj$theta, ni = tobj$bi, qi = tobj$bi2, 
        cgh = cgh, var = vartheta, vred = vred, y = y, hmax = hakt, 
        mae = mae, lseq = c(0, lseq[-steps]), call = args)
    class(z) <- "aws.gaussian"
  returns an object (list) of class \code{aws.gaussian} with components
  \item{theta }{Array of dimension \code{dim(y)} containing the smoothed data. }
  \item{ni }{Array of dimension \code{dim(y)[1:3]} containing the voxelwise sum of
  weights.}
  \item{qi }{Array of dimension \code{dim(y)[1:3]} containing the voxelwise sum of
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, qlambda = NULL, qtau = NULL, lkern = "Triangle", 
    aggkern = "Uniform", sigma2 = NULL, hinit = NULL, hincr = NULL, 
    hmax = NULL, lseq = NULL, heta = NULL, u = NULL, graph = FALSE, 
    demo = FALSE, wghts = NULL, spmin = 0, spmax = 5, scorr = 0, 
    vwghts = 1) 
{
    IQRdiff <- function(y) IQR(diff(y))/1.908
    updtheta <- function(zobj, tobj, cpar, aggkern) {
        heta <- cpar$heta
        tau1 <- cpar$tau1
        tau2 <- cpar$tau2
        kstar <- cpar$kstar
        hakt <- zobj$hakt
        tau <- 2 * (tau1 + tau2 * max(kstar - log(hakt), 0))
        hakt <- zobj$hakt
        bi <- zobj$bi
        bi2 <- zobj$bi2
        thetanew <- zobj$ai/bi
        theta <- tobj$theta
        thetanew[tobj$fix] <- theta[tobj$fix]
        if (hakt > heta) {
            eta <- switch(aggkern, Uniform = as.numeric(zobj$bi0/tau * 
                (thetanew - theta)^2 > 1), Triangle = pmin(1, 
                zobj$bi0/tau * (thetanew - theta)^2), as.numeric(zobj$bi0/tau * 
                (thetanew - theta)^2 > 1))
        }
        else {
            eta <- rep(0, length(theta))
        }
        eta[tobj$fix] <- 1
        theta <- (1 - eta) * thetanew + eta * theta
        eta <- eta[1:length(bi)]
        bi <- (1 - eta) * bi + eta * tobj$bi
        bi2 <- (1 - eta) * bi2 + eta * tobj$bi2
        list(theta = theta, bi = bi, bi2 = bi2, eta = eta, fix = (eta == 
            1))
    }
    args <- match.call()
    cat("\nvaws3D: start smoothing\n")
    d <- 3
    dy <- dim(y)
    n1 <- dy[1]
    n2 <- dy[2]
    n3 <- dy[3]
    n <- n1 * n2 * n3
    if (length(dy) == d) {
        dim(y) <- dy <- c(dy, 1)
    }
    else if (length(dy) != d + 1) {
        stop("y has to be 3 or 4 dimensional")
    }
    dv <- dim(y)[d + 1]
    if (length(vwghts) > dv) 
        vwghts <- vwghts[1:dv]
    dv0 <- length(vwghts)
    mae <- NULL
    if (is.null(heta)) 
        heta <- max(2, hinit + 0.5)
    heta <- max(heta, 2)
    if (is.null(qlambda)) 
        qlambda <- 0.985
    if (qlambda < 0.9) 
        warning("Inappropriate value of qlambda")
    if (qlambda >= 1) {
        if (is.null(qtau)) 
            qtau <- 0.4
        if (qtau == 1) 
            tau1 <- 1e+50
        else tau1 <- qchisq(qtau, 1)
        if (aggkern == "Triangle") 
            tau1 <- 2.5 * tau1
        tau2 <- tau1/2
    }
    else {
        if (is.null(qtau)) 
            qtau <- 0.95
        if (qtau >= 1) {
            tau1 <- 1e+50
        }
        else {
            tau1 <- qchisq(qtau, 1)
        }
        if (aggkern == "Triangle") 
            tau1 <- 2.5 * tau1
        tau2 <- tau1/2
    }
    if (qlambda < 1) {
        vwghts <- vwghts/max(vwghts)
        df <- sum(vwghts^2)^2/sum(vwghts^4)
        lambda <- qchisq(qlambda, df)
    }
    else {
        lambda <- 1e+50
    }
    cpar <- list(heta = heta, tau1 = tau1, tau2 = tau2)
    cpar$kstar <- log(5)
    lkern <- switch(lkern, Triangle = 2, Quadratic = 3, Cubic = 4, 
        Uniform = 1, Gaussian = 5, 2)
    if (is.null(hinit) || hinit < 1) 
        hinit <- 1
    if (is.null(hmax)) 
        hmax <- 5
    if (lkern == 5) {
        hmax <- hmax * 0.42445 * 4
        hinit <- min(hinit, hmax)
    }
    if (is.null(hincr) || hincr <= 1) 
        hincr <- 1.25
    hincr <- hincr^(1/d)
    h0 <- rep(0, length(scorr))
    if (max(scorr) > 0) {
        h0 <- numeric(length(scorr))
        for (i in 1:length(h0)) h0[i] <- get.bw.gauss(scorr[i], 
            interv = 2)
        if (length(h0) < d) 
            h0 <- rep(h0[1], d)
        cat("Corresponding bandwiths for specified correlation:", 
            h0, "\n")
    }
    if (is.null(sigma2)) {
        sigma2 <- IQRdiff(as.vector(y))^2
        if (scorr[1] > 0) 
            sigma2 <- sigma2 * Varcor.gauss(h0) * Spatialvar.gauss(h0 * 
                c(1, wghts), 1e-05, d)
        cat("Estimated variance: ", signif(sigma2, 4), "\n")
    }
    if (length(sigma2) == 1) {
        lambda <- lambda * sigma2 * 2
        cpar$tau1 <- cpar$tau1 * sigma2 * 2
        cpar$tau2 <- cpar$tau2 * sigma2 * 2
    }
    else {
        if (length(sigma2) != n) 
            stop("sigma2 does not have length 1 or same length as y")
        lambda <- lambda * 2
        cpar$tau1 <- cpar$tau1 * 2
        cpar$tau2 <- cpar$tau2 * 2
        sigma2 <- 1/sigma2
    }
    if (demo && !graph) 
        graph <- TRUE
    if (is.null(wghts)) 
        wghts <- c(1, 1, 1)
    hinit <- hinit/wghts[1]
    hmax <- hmax/wghts[1]
    wghts <- (wghts[2:3]/wghts[1])
    tobj <- list(bi = rep(1, n), bi2 = rep(1, n), theta = y, 
        fix = rep(FALSE, n))
    zobj <- list(ai = y, bi0 = rep(1, n))
    biold <- rep(1, n)
    if (length(sigma2) == n) 
        vred <- rep(1, n)
    steps <- as.integer(log(hmax/hinit)/log(hincr) + 1)
    if (is.null(lseq)) 
        lseq <- c(1.75, 1.35, 1.2, 1.2, 1.2, 1.2)
    if (length(lseq) < steps) 
        lseq <- c(lseq, rep(1, steps - length(lseq)))
    lseq <- lseq[1:steps]
    k <- 1
    hakt <- hinit
    hakt0 <- hinit
    lambda0 <- lambda
    if (hinit > 1) 
        lambda0 <- 1e+50
    progress <- 0
    step <- 0
    total <- (hincr^(d * ceiling(log(hmax/hinit)/log(hincr))) - 
        1)/(hincr^d - 1)
    while (hakt <= hmax) {
        dlw <- (2 * trunc(hakt/c(1, wghts)) + 1)[1:d]
        if (scorr[1] >= 0.1) 
            lambda0 <- lambda0 * Spatialvar.gauss(hakt0/0.42445/4 * 
                c(1, wghts), h0 * c(1, wghts), d)/Spatialvar.gauss(h0 * 
                c(1, wghts), 1e-05, d)/Spatialvar.gauss(hakt0/0.42445/4 * 
                c(1, wghts), 1e-05, d)
        hakt0 <- hakt
        if (length(sigma2) == n) {
            zobj <- .Fortran("chaws", as.double(y), as.logical(tobj$fix), 
                as.double(sigma2), as.integer(n1), as.integer(n2), 
                as.integer(n3), as.integer(dv), as.integer(dv0), 
                hakt = as.double(hakt), as.double(lambda0), as.double(tobj$theta), 
                bi = as.double(tobj$bi), bi2 = double(n), bi0 = as.double(zobj$bi0), 
                vred = double(n), ai = as.double(zobj$ai), as.integer(lkern), 
                as.double(spmin), as.double(spmax), double(prod(dlw)), 
                as.double(wghts), as.double(vwghts), double(dv), 
                double(dv0), double(dv0), PACKAGE = "fmri")[c("bi", 
                "bi0", "bi2", "vred", "ai", "hakt")]
            vred[!tobj$fix] <- zobj$vred[!tobj$fix]
        }
        else {
            zobj <- .Fortran("caws", as.double(y), as.logical(tobj$fix), 
                as.integer(n1), as.integer(n2), as.integer(n3), 
                as.integer(dv), as.integer(dv0), hakt = as.double(hakt), 
                as.double(lambda0), as.double(tobj$theta), bi = as.double(tobj$bi), 
                bi2 = double(n), bi0 = as.double(zobj$bi0), ai = as.double(zobj$ai), 
                as.integer(lkern), as.double(spmin), as.double(spmax), 
                double(prod(dlw)), as.double(wghts), as.double(vwghts), 
                double(dv), double(dv0), double(dv0), PACKAGE = "fmri")[c("bi", 
                "bi0", "bi2", "ai", "hakt")]
        }
        gc()
        dim(zobj$ai) <- dy
        if (hakt > n1/2) 
            zobj$bi0 <- hincr^d * biold
        biold <- zobj$bi0
        tobj <- updtheta(zobj, tobj, cpar, aggkern)
        gc()
        dim(tobj$theta) <- dy
        dim(tobj$bi) <- dy[-4]
        dim(tobj$eta) <- dy[-4]
        if (graph) {
            par(mfrow = c(2, 2), mar = c(1, 1, 3, 0.25), mgp = c(2, 
                1, 0))
            image(y[, , n3\%/\%2 + 1, 1], col = gray((0:255)/255), 
                xaxt = "n", yaxt = "n")
            title(paste("Observed Image  min=", signif(min(y), 
                3), " max=", signif(max(y), 3)))
            image(tobj$theta[, , n3\%/\%2 + 1, 1], col = gray((0:255)/255), 
                xaxt = "n", yaxt = "n")
            title(paste("Reconstruction  h=", signif(hakt, 3), 
                " min=", signif(min(tobj$theta), 3), " max=", 
                signif(max(tobj$theta), 3)))
            image(tobj$bi[, , n3\%/\%2 + 1], col = gray((0:255)/255), 
                xaxt = "n", yaxt = "n")
            title(paste("Sum of weights: min=", signif(min(tobj$bi), 
                3), " mean=", signif(mean(tobj$bi), 3), " max=", 
                signif(max(tobj$bi), 3)))
            image(tobj$eta[, , n3\%/\%2 + 1], col = gray((0:255)/255), 
                xaxt = "n", yaxt = "n", zlim = c(0, 1))
            title("eta")
        }
        if (!is.null(u)) {
            cat("bandwidth: ", signif(hakt, 3), "eta==1", sum(tobj$eta == 
                1), "   MSE: ", signif(mean((tobj$theta - u)^2), 
                3), "   MAE: ", signif(mean(abs(tobj$theta - 
                u)), 3), " mean(bi)=", signif(mean(tobj$bi), 
                3), "\n")
            mae <- c(mae, signif(mean(abs(tobj$theta - u)), 3))
        }
        else if (total != 0) {
            progress <- progress + hincr^(d * step)
            step <- step + 1
            cat(signif(progress/total, 2) * 100, "\% . ", sep = "")
        }
        if (demo) 
            readline("Press return")
        hakt <- hakt * hincr
        x <- 1.25^(k - 1)
        scorrfactor <- x/(3^d * prod(scorr) * prod(h0) + x)
        lambda0 <- lambda * lseq[k] * scorrfactor
        k <- k + 1
        gc()
    }
    if (length(sigma2) == n) {
        vartheta <- tobj$bi2/tobj$bi^2
    }
    else {
        vartheta <- sigma2 * tobj$bi2/tobj$bi^2
        vred <- tobj$bi2/tobj$bi^2
    }
    hakt <- hakt/hincr
    cgh <- Spatialvar.gauss(hakt/0.42445/4 * c(1, wghts), h0 * 
        c(1, wghts) + 1e-05, d)/Spatialvar.gauss(hakt/0.42445/4 * 
        c(1, wghts), 1e-05, d)/Spatialvar.gauss(h0 * c(1, wghts), 
        1e-05, d)
    vartheta <- vartheta * cgh
    z <- list(theta = tobj$theta, ni = tobj$bi, qi = tobj$bi2, 
        cgh = cgh, var = vartheta, vred = vred, y = y, hmax = hakt, 
        mae = mae, lseq = c(0, lseq[-steps]), call = args)
    class(z) <- "aws.gaussian"
    cat("\nvaws3D: finished smoothing\n")
    z
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
