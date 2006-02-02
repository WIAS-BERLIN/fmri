\name{perform.aws}
\alias{perform.aws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
perform.aws(beta, variance, hmax = 4, hinit = 1, weights = c(1, 1, 1), qlambda = 1, lkern = "Gaussian", scorr = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{ ~~Describe \code{beta} here~~ }
  \item{variance}{ ~~Describe \code{variance} here~~ }
  \item{hmax}{ ~~Describe \code{hmax} here~~ }
  \item{hinit}{ ~~Describe \code{hinit} here~~ }
  \item{weights}{ ~~Describe \code{weights} here~~ }
  \item{qlambda}{ ~~Describe \code{qlambda} here~~ }
  \item{lkern}{ ~~Describe \code{lkern} here~~ }
  \item{scorr}{ ~~Describe \code{scorr} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
function (beta, variance, hmax = 4, hinit = 1, weights = c(1, 
    1, 1), qlambda = 1, lkern = "Gaussian", scorr = 0, ...) 
{
    variance[variance < quantile(variance, 0.25)] <- quantile(variance, 
        0.25)
    ttthat <- vaws3D(beta, sigma2 = variance, hmax = hmax, hinit = hinit, 
        qlambda = qlambda, qtau = 1, heta = 1000, wghts = weights, 
        lkern = lkern, scorr = scorr, ...)
    if (dim(ttthat$theta)[4] == 1) {
        list(hat = ttthat$theta[, , , 1], var = ttthat$var, ni = ttthat$ni, 
            hmax = ttthat$hmax, lseq = ttthat$lseq, y = ttthat$y, 
            mae = ttthat$mae, vred = ttthat$vred)
    }
    else {
        list(hat = ttthat$theta, var = ttthat$var, ni = ttthat$ni, 
            hmax = ttthat$hmax, lseq = ttthat$lseq, y = ttthat$y, 
            mae = ttthat$mae, vred = ttthat$vred)
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
