\name{calculate.lm}
\alias{calculate.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculate.lm(ttt, z, actype = "smooth", hmax = 3.52, vtype = "var", step = 0.01, contrast = c(1), vvector = c(1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttt}{ ~~Describe \code{ttt} here~~ }
  \item{z}{ ~~Describe \code{z} here~~ }
  \item{actype}{ ~~Describe \code{actype} here~~ }
  \item{hmax}{ ~~Describe \code{hmax} here~~ }
  \item{vtype}{ ~~Describe \code{vtype} here~~ }
  \item{step}{ ~~Describe \code{step} here~~ }
  \item{contrast}{ ~~Describe \code{contrast} here~~ }
  \item{vvector}{ ~~Describe \code{vvector} here~~ }
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
function (ttt, z, actype = "smooth", hmax = 3.52, vtype = "var", 
    step = 0.01, contrast = c(1), vvector = c(1)) 
{
    create.arcorrection <- function(scans, rho = 0) {
        rho0 <- 1/sqrt(1 - rho^2)
        a <- numeric(scans * scans)
        a[1] <- 1
        ind <- (2:scans) * (scans + 1) - 2 * scans
        a[ind] <- -rho * rho0
        a[ind + scans] <- rho0
        dim(a) <- c(scans, scans)
        a
    }
    cat("calculate.lm: entering function with:", actype, hmax, 
        vtype, "\n")
    if (length(contrast) <= dim(z)[2]) 
        contrast <- c(contrast, rep(0, dim(z)[2] - length(contrast)))
    length(contrast) <- dim(z)[2]
    if ((length(vvector) < dim(z)[2]) && (length(vvector) > 1)) 
        vvector <- c(vvector, rep(0, dim(z)[2] - length(vvector)))
    svdresult <- svd(z)
    u <- svdresult$u
    v <- svdresult$v
    vt <- t(v)
    lambda1 <- diag(1/svdresult$d)
    lambda2 <- diag(1/svdresult$d^2)
    xtx <- v \%*\% lambda2 \%*\% vt
    dy <- dim(ttt)
    voxelcount <- prod(dy[1:3])
    dim(ttt) <- c(prod(dy[1:3]), dy[4])
    arfactor <- rep(0, length = prod(dy[1:3]))
    variance <- rep(0, length = prod(dy[1:3]))
    if (length(vvector) > 1) 
        variancem <- array(0, c(dy[1:3], length(vector)^2))
    R <- diag(1, dy[4]) - u \%*\% t(u)
    m00 <- dy[4] - dim(z)[2]
    m01 <- 0
    for (k in 1:dy[4]) {
        m01 <- m01 + sum(R[k, -1] * R[k, -dy[4]])
    }
    m10 <- m01
    m11 <- 0
    for (k in 1:(dy[4] - 1)) {
        m11 <- m11 + sum(R[k + 1, -dy[4]] * R[k, -1] + R[k + 
            1, -1] * R[k, -dy[4]])
    }
    Minv <- matrix(c(m11, -m01, -m10, m00), 2, 2)/(m00 * m11 - 
        m01 * m10)
    beta <- ttt \%*\% u \%*\% lambda1 \%*\% vt
    residuals <- ttt - beta \%*\% t(z)
    if ((actype == "smooth") || (actype == "accalc") || (actype == 
        "ac")) {
        progress = 0
        cat("calculate.lm: calculating AR(1) model\n")
        for (i in 1:voxelcount) {
            if (i > progress/100 * voxelcount) {
                cat(progress, "\% . ", sep = "")
                progress = progress + 10
            }
            a0 <- residuals[i, ] \%*\% residuals[i, ]
            a1 <- residuals[i, -1] \%*\% residuals[i, -dim(z)[1]]
            an <- Minv \%*\% c(a0, a1)
            if (an[1] != 0) {
                arfactor[i] <- an[2]/an[1]
            }
            else {
                arfactor[i] <- 0
            }
        }
        cat("\n")
        cat("calculate.lm: finished\n")
        if (actype == "smooth") {
            cat("calculate.lm: smoothing with (hmax):", hmax, 
                "\n")
            dim(arfactor) <- dy[1:3]
            hinit <- 1
            arfactor <- gkernsm(arfactor, rep(hmax, 3) * 0.42445)$gkernsm
            dim(arfactor) <- voxelcount
            cat("calculate.lm: finished\n")
        }
        if ((actype == "smooth") || (actype == "accalc")) {
            progress = 0
            cat("calculate.lm: re-calculating linear model with prewithened data\n")
            variancepart <- rep(0, length = prod(dy[1:3]))
            if (length(vvector) > 1) 
                variancepartm <- array(0, c(sum(as.logical(vvector))^2, 
                  prod(dy[1:3])))
            arlist <- seq(range(arfactor)[1] - step/2, range(arfactor)[2] + 
                step/2, step)
            for (i in 1:(length(arlist) - 1)) {
                if (i > progress/100 * length(arlist)) {
                  cat(progress, "\% . ", sep = "")
                  progress = progress + 10
                }
                indar <- as.logical((arfactor > arlist[i]) * 
                  (arfactor <= arlist[i + 1]))
                if (sum(indar) > 0) {
                  a <- create.arcorrection(dy[4], mean(arlist[i:(i + 
                    1)]))
                  zprime <- a \%*\% z
                  svdresult <- svd(zprime)
                  v <- svdresult$v
                  vt <- t(v)
                  xtx <- v \%*\% diag(1/svdresult$d^2) \%*\% vt
                  tttprime <- ttt[indar, ] \%*\% t(a)
                  beta[indar, ] <- tttprime \%*\% svdresult$u \%*\% 
                    diag(1/svdresult$d) \%*\% vt
                  residuals[indar, ] <- tttprime - beta[indar, 
                    ] \%*\% t(zprime)
                  variancepart[indar] <- t(contrast) \%*\% xtx \%*\% 
                    contrast
                  if (length(vvector) > 1) 
                    variancepartm[, indar] <- as.vector(xtx[as.logical(vvector), 
                      as.logical(vvector)])
                }
            }
            b <- rep(1/dy[4], length = dy[4])
            variance <- ((residuals^2 \%*\% b) * dim(z)[1]/(dim(z)[1] - 
                dim(z)[2])) * variancepart
            if (length(vvector) > 1) 
                variancem <- as.vector((residuals^2 \%*\% b) * 
                  dim(z)[1]/(dim(z)[1] - dim(z)[2])) * t(variancepartm)
            cat("\n")
            cat("calculate.lm: finished\n")
        }
        else {
            b <- rep(1/dy[4], length = dy[4])
            cxtx <- t(contrast) \%*\% xtx \%*\% contrast
            if (length(vvector) > 1) 
                variancem <- ((residuals^2 \%*\% b) * dim(z)[1]/(dim(z)[1] - 
                  dim(z)[2])) \%*\% as.vector(xtx[as.logical(vvector), 
                  as.logical(vvector)])
            variance <- ((residuals^2 \%*\% b) * dy[4]/(dy[4] - 
                dim(z)[2])) \%*\% cxtx
        }
    }
    else {
        b <- rep(1/dy[4], length = dy[4])
        cxtx <- t(contrast) \%*\% xtx \%*\% contrast
        if (length(vvector) > 1) 
            variancem <- ((residuals^2 \%*\% b) * dim(z)[1]/(dim(z)[1] - 
                dim(z)[2])) \%*\% as.vector(xtx[as.logical(vvector), 
                as.logical(vvector)])
        variance <- ((residuals^2 \%*\% b) * dy[4]/(dy[4] - dim(z)[2])) \%*\% 
            cxtx
    }
    cbeta <- beta \%*\% contrast
    dim(beta) <- c(dy[1:3], dim(z)[2])
    dim(cbeta) <- dy[1:3]
    dim(variance) <- dy[1:3]
    if (length(vvector) > 1) 
        dim(variancem) <- c(dy[1:3], sum(as.logical(vvector)), 
            sum(as.logical(vvector)))
    dim(arfactor) <- dy[1:3]
    dim(residuals) <- dy
    cat("calculate.lm: exiting function\n")
    if (length(vvector) > 1) {
        list(beta = beta, cbeta = cbeta, var = variance, varm = variancem, 
            res = residuals, arfactor = arfactor)
    }
    else {
        list(beta = beta, cbeta = cbeta, var = variance, res = residuals, 
            arfactor = arfactor)
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
