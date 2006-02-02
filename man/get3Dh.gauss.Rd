\name{get3Dh.gauss}
\alias{get3Dh.gauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get3Dh.gauss(vred, h0, vwghts, step = 1.002, interv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vred}{ ~~Describe \code{vred} here~~ }
  \item{h0}{ ~~Describe \code{h0} here~~ }
  \item{vwghts}{ ~~Describe \code{vwghts} here~~ }
  \item{step}{ ~~Describe \code{step} here~~ }
  \item{interv}{ ~~Describe \code{interv} here~~ }
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
function (vred, h0, vwghts, step = 1.002, interv = 1) 
{
    h0 <- pmax(h0, 1e-05)
    n <- length(vred)
    vred1 <- vred
    h <- 0.5/vwghts
    fixed <- rep(FALSE, length(vred))
    while (any(!fixed)) {
        ind <- (1:n)[!fixed][vred[!fixed] >= Spatialvar.gauss(h, 
            1e-05, 3, interv)]
        vred1[ind] <- Spatialvar.gauss(h, h0, 3, interv)
        fixed[ind] <- TRUE
        h <- h * step
    }
    hvred <- matrix(0, 3, n)
    hh <- 0.01/vwghts
    h <- h0
    fixed <- rep(FALSE, length(vred))
    while (any(!fixed)) {
        ind <- (1:n)[!fixed][vred1[!fixed] >= Spatialvar.gauss(h, 
            1e-05, 3, interv)]
        hvred[, ind] <- h
        fixed[ind] <- TRUE
        hh <- hh * step
        h <- sqrt(h0^2 + hh^2)
    }
    t(hvred)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
