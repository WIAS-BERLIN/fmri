\name{vaws3D}
\alias{vaws3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vaws3D(y, qlambda = NULL, qtau = NULL, lkern = "Triangle", aggkern = "Uniform", sigma2 = NULL, hinit = NULL, hincr = NULL, hmax = NULL, lseq = NULL, heta = NULL, u = NULL, graph = FALSE, demo = FALSE, wghts = NULL, spmin = 0, spmax = 5, scorr = 0, vwghts = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{qlambda}{ ~~Describe \code{qlambda} here~~ }
  \item{qtau}{ ~~Describe \code{qtau} here~~ }
  \item{lkern}{ ~~Describe \code{lkern} here~~ }
  \item{aggkern}{ ~~Describe \code{aggkern} here~~ }
  \item{sigma2}{ ~~Describe \code{sigma2} here~~ }
  \item{hinit}{ ~~Describe \code{hinit} here~~ }
  \item{hincr}{ ~~Describe \code{hincr} here~~ }
  \item{hmax}{ ~~Describe \code{hmax} here~~ }
  \item{lseq}{ ~~Describe \code{lseq} here~~ }
  \item{heta}{ ~~Describe \code{heta} here~~ }
  \item{u}{ ~~Describe \code{u} here~~ }
  \item{graph}{ ~~Describe \code{graph} here~~ }
  \item{demo}{ ~~Describe \code{demo} here~~ }
  \item{wghts}{ ~~Describe \code{wghts} here~~ }
  \item{spmin}{ ~~Describe \code{spmin} here~~ }
  \item{spmax}{ ~~Describe \code{spmax} here~~ }
  \item{scorr}{ ~~Describe \code{scorr} here~~ }
  \item{vwghts}{ ~~Describe \code{vwghts} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
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
