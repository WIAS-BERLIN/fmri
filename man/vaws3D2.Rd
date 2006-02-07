\name{vaws3D2}
\alias{vaws3D2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Three dimensional adaptive smoothing of vector valued data}
\description{Performs three dimensional adaptive smoothing of vector valued data
using the Propagation-Separation approach from Polzehl \& Spokoiny (2006).  }
\usage{
vaws3D2(y, qlambda = NULL, qtau = NULL, lkern = "Triangle", skern = "Exp", aggkern = "Uniform", sigma2 = NULL, hinit = NULL, hincr = NULL, hmax = NULL, lseq = NULL, heta = NULL, u = NULL, graph = FALSE, demo = FALSE, wghts = NULL, spmin = 0, spmax = 1.1, scorr = c(0, 0, 0), vwghts = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Data array of dimension \code{c(dx,dy,dz,dv)}. The
  first three values specify the size of the 3D data volume. \code{dv}
  specifies the length of the data vector in each voxel.  }
  \item{qlambda}{\code{qlambda} is the main smoothing parameter. It determines the scaling parameter \code{\lambda} in the statistical penalty
  (see description of the algorithm in the references) as \code{qchisq(qlambda,sum(vwghts^2)^2/sum(vwghts^4))}. \code{qlambda=1} disables adaptation, i.e. the 
  resulting estimate is a kernel estimate using the largest inspected bandwidth less or equal to \code{hmax}.
  The default value for \code{qlambda} is selected as the smallest value that fulfils a propagation condition.
  This means that in a parametric, in this case constant, model the resulting estimate is, for a large \code{hmax} almost parametric. 
   Larger values of 
  \code{qlambda} lead to less sensitivity to structural differences, smaller values may lead to a random structure (segmentation) of the resulting estimate.}
  \item{qtau}{ not used anymore to be removed }
  \item{lkern}{ \code{lkern} specifies the location kernel. Defaults to "Triangle", other choices are "Gaussian", "Quadratic", "Cubic" and "Uniform".
    Note that the location kernel is applied to \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression. }
  \item{skern}{\code{lkern} specifies the kernel for the statistikal penalty. Defaults to "Exp", the alternative is "Triangle".
    \code{lkern="Triangle"} allows for much faster computation (saves up to 50%) but needs some additional testing. }
  \item{aggkern}{\code{aggkern} specifies the kernel for the statistical panalty in stagewise aggregation. Defaults to "Uniform", the alternative choice is
   "Triangle"  }
  \item{sigma2}{\code{sigma2} allows to specify the variance \eqn{\sigma^2_i} of the first component of the data vector in each voxel \eqn{i}. 
   Defaults to \code{NULL}. In this case a homoskedastic variance estimate is generated. If \code{length(sigma2)==1} a homoskedastic
   variance is assumed.  }
  \item{hinit}{\code{hinit} specifies the initial bandwidth. Defaults to \code{hinit=1} (\code{hinit=1/4} in case of \code{lkern="Gaussian"}).   }
  \item{hincr}{\code{hincr} specifies the factor used to increase the size of local neigborhoods after each iteration. The bandwidth is increased by
   a factor \code{hincr^(1/3)}. Defaults to \code{hincr=1.25}   }
  \item{hmax}{\code{hmax} specifies the maximal bandwidth. Defaults to \code{hmax= 5}}
  \item{lseq}{ \code{lseq} allows to increase the value of the scaling parameter \code{\lambda} for the first \code{length(lseq)} iteration steps by
  the factor specified in \code{lseq}. Defaults to NULL. In this case a default sequence is used that fulfils the propagation condition.}
  \item{heta}{ not used anymore to be removed }
  \item{u}{\code{u} allows to specify the true parameter. This is only used to test the algorithm and to select the smoothing parameters 
    \code{qlambda} and \code{qtau} by a propagation condition. If \code{u} is specified MSE and MAE of the estimates are 
    printed for each iteration step.  }
  \item{graph}{If  \code{graph=TRUE} intermediate results are illustrated after each iteration step. Defaults to \code{graph=FALSE}.  }
  \item{demo}{If \code{demo=TRUE} the function pauses after each iteration. Defaults to \code{demo=FALSE}. }
  \item{wghts}{\code{wghts} specifies the  diagonal elements of a weight matrix to adjust for different distances between grid-points
  in different coordinate directions, i.e. allows to define a more appropriate metric in the design space.  }
  \item{spmin}{ \code{spmin} specifies parameter value in the kernel \eqn{K_s(x)=\min(0,e^{-spmax/(spmax-spmin)(x-spmin)}I_{x<spmax}}. Defaults
  to \code{spmin=0}. Not used if \code{skern="Triangle"}.}
  \item{spmax}{\code{spmax} specifies the cut point in the kernel \eqn{K_s(x)=\min(0,e^{-spmax/(spmax-spmin)(x-spmin)}I_{x<spmax}} }. Defaults
  to \code{spmax=5}. Not used if \code{skern="Triangle"}.}
  \item{scorr}{\code{scorr} specifies the spatial correlation in the three coordinate directions. Defaults to \code{scorr=c(0,0,0)}. }
  \item{vwghts}{\code{vwghts} specifies weights for the components of the data vector. These weights are used in the evaluation of the statistical
  penalty \eqn{s_ij=N_i/\lambda/\sigma^2_i sum_k vwghts[k] (theta_{ik}-\theta_{jk})^2}{} where \eqn{k} denotes the \eqn{k}-th component of the 
  vector. The parameter can be used to adjust for different variability of the components of the data vector. Defaults to \code{vwghts=rep(1,dv)}.  }
}
\details{
 This function differs from \code{\link{vaws3D}} in three points:\cr
 1) it computes  estimates \code{vartheta} and \code{vred} (variance reduction)
 that are more accurate.\cr
 2) it does not implement stagewise aggregation, that is, corresponds to the
 case \code{qtau=1}.\cr
 3) the information needed for \code{vartheta} and \code{vred} is only
 evaluated at after the estimates are fixed.\cr
 
 The function implements the propagation separation approach to 
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
   \eqn{\bf{E} |\hat{\theta}^k - \theta| \le (1+\alpha) \bf{E} |\tilde{\theta}^k - \theta|}
   where \code{\hat{\theta}^k} is the aws-estimate in step \code{k} and \eqn{\tilde{\theta}^k}
   is corresponding nonadaptive estimate using the same bandwidth (\code{qlambda=1}).
   Default values are selected to fulfil this condition for \eqn{\alpha=0.1}.
   
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
  \item{var}{The estimated variances of \code{theta[,,,1]}.}
  \item{vred}{The amount of variance reduction achieved in each voxel.}
  \item{y}{The data.}
  \item{qi}{Array of dimension \code{dim(y)[1:3]} containing the voxelwise sum of
  squared weights.}
  \item{cgh}{The value \code{cgh}.}
  \item{hmax}{The largest bandwidth used.}
  \item{mae}{Mean absolute error if \code{!is.null(u)}.}
  \item{lseq}{The sequence specified in \code{lseq}.}
  \item{call}{The arguments used in the call of function vaws3D}  
  \item{ng}{Array of dimension \code{dim(y)[1:3]} containing the voxelwise sums \eqn{N_g}
  \item{qg}{Array of dimension \code{dim(y)[1:3]} containing the voxelwise sums \eqn{Q_g}
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, qlambda = NULL, qtau = NULL, lkern = "Triangle", 
    skern = "Exp", aggkern = "Uniform", sigma2 = NULL, hinit = NULL, 
    hincr = NULL, hmax = NULL, lseq = NULL, heta = NULL, u = NULL, 
    graph = FALSE, demo = FALSE, wghts = NULL, spmin = 0, spmax = 1.1, 
    scorr = c(0, 0, 0), vwghts = 1) 
{
    IQRdiff <- function(y) IQR(diff(y))/1.908
    args <- match.call()
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
    lkern <- switch(lkern, Triangle = 2, Quadratic = 3, Cubic = 4, 
        Uniform = 1, Gaussian = 5, 2)
    skern <- switch(skern, Triangle = 1, Exp = 2, 2)
    if (is.null(qlambda)) 
        qlambda <- 0.985
    if (qlambda < 0.9) 
        warning("Inappropriate value of qlambda")
    if (qlambda < 1) {
        vwghts <- vwghts/max(vwghts)
        df <- sum(vwghts^2)^2/sum(vwghts^4)
        lambda <- qchisq(qlambda, df)
    }
    else {
        lambda <- 1e+50
    }
    if (skern == 1) {
        lambda <- 4/3 * lambda
        spmin <- 0
        spmax <- 1
    }
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
    if (length(sigma2) == 1) 
        sigma2 <- array(sigma2, dy[1:3])
    if (length(sigma2) != n) 
        stop("sigma2 does not have length 1 or same length as y")
    lambda <- lambda * 2
    sigma2 <- 1/sigma2
    if (demo && !graph) 
        graph <- TRUE
    if (is.null(wghts)) 
        wghts <- c(1, 1, 1)
    hinit <- hinit/wghts[1]
    hmax <- hmax/wghts[1]
    wghts <- (wghts[2:3]/wghts[1])
    tobj <- list(ai = y, bi = rep(1, n))
    theta <- y
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
    while (hakt <= hmax) {
        dlw <- (2 * trunc(hakt/c(1, wghts)) + 1)[1:d]
        if (scorr[1] >= 0.1) 
            lambda0 <- lambda0 * Spatialvar.gauss(hakt0/0.42445/4 * 
                c(1, wghts), h0 * c(1, wghts), d)/Spatialvar.gauss(h0 * 
                c(1, wghts), 1e-05, d)/Spatialvar.gauss(hakt0/0.42445/4 * 
                c(1, wghts), 1e-05, d)
        hakt0 <- hakt
        theta0 <- theta
        bi0 <- tobj$bi
        tobj <- .Fortran("chaws2", as.double(y), as.double(sigma2), 
            as.integer(n1), as.integer(n2), as.integer(n3), as.integer(dv), 
            as.integer(dv0), hakt = as.double(hakt), as.double(lambda0), 
            as.double(theta0), bi = as.double(bi0), ai = as.double(tobj$ai), 
            as.integer(lkern), as.integer(skern), as.double(spmin), 
            as.double(spmax), double(prod(dlw)), as.double(wghts), 
            as.double(vwghts), double(dv), double(dv0), double(dv0), 
            PACKAGE = "fmri", DUP = TRUE)[c("bi", "ai", "hakt")]
        gc()
        dim(tobj$ai) <- dy
        gc()
        theta <- array(tobj$ai/tobj$bi, dy)
        dim(tobj$bi) <- dy[-4]
        if (graph) {
            par(mfrow = c(1, 3), mar = c(1, 1, 3, 0.25), mgp = c(2, 
                1, 0))
            image(y[, , n3\%/\%2 + 1, 1], col = gray((0:255)/255), 
                xaxt = "n", yaxt = "n")
            title(paste("Observed Image  min=", signif(min(y), 
                3), " max=", signif(max(y), 3)))
            image(theta[, , n3\%/\%2 + 1, 1], col = gray((0:255)/255), 
                xaxt = "n", yaxt = "n")
            title(paste("Reconstruction  h=", signif(hakt, 3), 
                " min=", signif(min(theta), 3), " max=", signif(max(theta), 
                  3)))
            image(tobj$bi[, , n3\%/\%2 + 1], col = gray((0:255)/255), 
                xaxt = "n", yaxt = "n")
            title(paste("Sum of weights: min=", signif(min(tobj$bi), 
                3), " mean=", signif(mean(tobj$bi), 3), " max=", 
                signif(max(tobj$bi), 3)))
        }
        if (!is.null(u)) {
            cat("bandwidth: ", signif(hakt, 3), "eta==1", sum(tobj$eta == 
                1), "   MSE: ", signif(mean((theta - u)^2), 3), 
                "   MAE: ", signif(mean(abs(theta - u)), 3), 
                " mean(bi)=", signif(mean(tobj$bi), 3), "\n")
            mae <- c(mae, signif(mean(abs(theta - u)), 3))
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
    g <- trunc(h0/c(1, wghts)/2.3548 * 4)
    if (h0[1] > 0) 
        gw1 <- dnorm(-(g[1]):g[1], 0, h0[1]/2.3548)/dnorm(0, 
            0, h0[1]/2.3548)
    else gw1 <- 1
    if (h0[2] > 0) 
        gw2 <- dnorm(-(g[2]):g[2], 0, h0[2]/2.3548/wghts[1])/dnorm(0, 
            0, h0[2]/2.3548/wghts[1])
    else gw2 <- 1
    if (h0[3] > 0) 
        gw3 <- dnorm(-(g[3]):g[3], 0, h0[3]/2.3548/wghts[2])/dnorm(0, 
            0, h0[3]/2.3548/wghts[2])
    else gw3 <- 1
    gwght <- outer(gw1, outer(gw2, gw3, "*"), "*")
    gwght <- gwght/sum(gwght)
    dgw <- dim(gwght)
    z <- .Fortran("chawsvr", as.double(y), as.double(sigma2), 
        as.integer(n1), as.integer(n2), as.integer(n3), as.integer(dv), 
        as.integer(dv0), as.double(tobj$hakt), as.double(lambda0), 
        as.double(theta0), as.double(bi0), var = double(n), vred = double(n), 
        as.integer(lkern), as.integer(skern), as.double(spmin), 
        as.double(spmax), double(prod(dlw)), as.double(gwght), 
        double(n), as.integer(dgw), as.double(wghts), as.double(vwghts), 
        double(dv0), double(dv0), PACKAGE = "fmri", DUP = FALSE)[c("vred", 
        "var")]
    dim(z$vred) <- dim(z$var) <- dim(sigma2)
    nqg <- .Fortran("nqg", as.double(gwght), as.double(gwght^2), 
        as.integer(dgw[1]), as.integer(dgw[2]), as.integer(dgw[3]), 
        as.integer(n1), as.integer(n2), as.integer(n3), qg = double(n), 
        ng = double(n), PACKAGE = "fmri", DUP = FALSE)[c("qg", 
        "ng")]
    qg <- array(nqg$qg, dim(sigma2))
    ng <- array(nqg$ng, dim(sigma2))
    dim(tobj$bi) <- dim(sigma2)
    vartheta <- z$var/tobj$bi^2/qg
    vred <- z$vred/tobj$bi^2/ng^2
    z <- list(theta = theta, ni = tobj$bi, var = vartheta, vred = vred, 
        y = y, hmax = tobj$hakt, mae = mae, lseq = c(0, lseq[-steps]), 
        call = args, ng = ng, qg = qg, gw = gwght)
    class(z) <- "aws.gaussian"
    z
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
