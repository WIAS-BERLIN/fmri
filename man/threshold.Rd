\name{threshold}
\alias{threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
threshold(p, i, j, k, rx, ry, rz, type = "norm", df = 4, dim = 3, step = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ ~~Describe \code{p} here~~ }
  \item{i}{ ~~Describe \code{i} here~~ }
  \item{j}{ ~~Describe \code{j} here~~ }
  \item{k}{ ~~Describe \code{k} here~~ }
  \item{rx}{ ~~Describe \code{rx} here~~ }
  \item{ry}{ ~~Describe \code{ry} here~~ }
  \item{rz}{ ~~Describe \code{rz} here~~ }
  \item{type}{ ~~Describe \code{type} here~~ }
  \item{df}{ ~~Describe \code{df} here~~ }
  \item{dim}{ ~~Describe \code{dim} here~~ }
  \item{step}{ ~~Describe \code{step} here~~ }
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
function (p, i, j, k, rx, ry, rz, type = "norm", df = 4, dim = 3, 
    step = 0.001) 
{
    n <- length(rx)
    thr <- numeric(n)
    fixed <- logical(n)
    x <- 3
    pxyz <- rx * ry * rz
    ind <- (1:length(pxyz))[pxyz == min(pxyz)][1]
    pv <- 1
    while (pv > p) {
        pv <- pvalue(x, i, j, k, rx[ind], ry[ind], rz[ind], type, 
            df, dim)
        x <- x + 5 * step
    }
    x <- x - 5 * step
    while (any(!fixed)) {
        pv <- pvalue(x, i, j, k, rx, ry, rz, type, df, dim)
        ind <- (1:n)[!fixed][pv[!fixed] < p]
        thr[ind] <- x
        fixed[ind] <- TRUE
        x <- x + step
    }
    thr
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
