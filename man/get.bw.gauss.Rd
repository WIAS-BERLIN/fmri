\name{get.bw.gauss}
\alias{get.bw.gauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.bw.gauss(corr, step = 1.001, interv = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corr}{ ~~Describe \code{corr} here~~ }
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
function (corr, step = 1.001, interv = 2) 
{
    get.corr.gauss <- function(h, interv = 1) {
        h <- h/2.3548 * interv
        ih <- trunc(4 * h + 2 * interv - 1)
        dx <- 2 * ih + 1
        penl <- dnorm(((-ih):ih)/h)
        sum(penl[-(1:interv)] * penl[-((dx - interv + 1):dx)])/sum(penl^2)
    }
    if (corr < 0.1) {
        h <- 0
    }
    else {
        h <- 0.5
        z <- 0
        while (z < corr) {
            h <- h * step
            z <- get.corr.gauss(h, interv)
        }
        h <- h/step
    }
    h
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
