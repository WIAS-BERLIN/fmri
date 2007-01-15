\name{ngca}
\alias{ngca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Non-Gaussian Component Analysis }
\description{
  The function performs Non-Gaussian Component Analysis as described in Blanchard et.al. (2006).
}
\usage{
ngca(x, L = 1000, T = 10, m = 3, eps = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Observation matrix (dimension Nxd) }
  \item{L}{ Number basis functions in each of four classes. }
  \item{T}{ Number of Fast ICA iterations  }
  \item{m}{ Number of non-Gaussian components. }
  \item{eps}{ Threshold (defaults to 1.5)}
}
\details{
  The function performs Non-Gaussian Component Analysis as described in Blanchard et.al. (2006).
  The procedure  classes of basis functions 
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
function (x, L = 1000, T = 10, m = 3, eps = 1.5) 
{
    xdim <- dim(x)
    d <- xdim[2]
    n <- xdim[1]
    xmean <- apply(x, 2, mean)
    xvar <- var(x)
    y <- sweep(x, 2, xmean)
    z <- svd(xvar)
    y <- y \%*\% z$u \%*\% diag(z$d^(-0.5))
    y <- t(y)
    s <- matrix(0, L, 4)
    s[, 1] <- seq(0.5, 5, length = L)
    s[, 2] <- seq(5/L, 5, length = L)
    s[, 3] <- seq(4/L, 4, length = L)
    s[, 4] <- seq(0, 4, length = L)
    omega <- matrix(rnorm(4 * L * d), d, L * 4)
    omega <- sweep(omega, 2, sqrt(apply(omega^2, 2, sum)), "/")
    fz <- .Fortran("fastica", as.double(y), as.double(omega), 
        as.integer(d), as.integer(n), as.integer(L), as.integer(T), 
        double(d), v = double(d * L * 4), normv = double(L * 
            4), as.double(s), DUP = FALSE, PACKAGE = "fmri")[c("v", 
        "normv")]
    dim(fz$v) <- c(d, 4 * L)
    fz$v <- t(fz$v[, fz$normv > eps])
    jhat <- prcomp(fz$v)
    ihat <- z$u \%*\% jhat$rotation[, 1:m]
    list(ihat = ihat, sdev = jhat$sdev[1:m])
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
