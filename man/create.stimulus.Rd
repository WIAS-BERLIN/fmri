\name{create.stimulus}
\alias{create.stimulus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create.stimulus(scans = 1, onsets = c(1), length = 1, rt = 3, mean = TRUE, a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scans}{ ~~Describe \code{scans} here~~ }
  \item{onsets}{ ~~Describe \code{onsets} here~~ }
  \item{length}{ ~~Describe \code{length} here~~ }
  \item{rt}{ ~~Describe \code{rt} here~~ }
  \item{mean}{ ~~Describe \code{mean} here~~ }
  \item{a1}{ ~~Describe \code{a1} here~~ }
  \item{a2}{ ~~Describe \code{a2} here~~ }
  \item{b1}{ ~~Describe \code{b1} here~~ }
  \item{b2}{ ~~Describe \code{b2} here~~ }
  \item{cc}{ ~~Describe \code{cc} here~~ }
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
function (scans = 1, onsets = c(1), length = 1, rt = 3, mean = TRUE, 
    a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35) 
{
    numberofonsets <- length(onsets)
    stimulus <- rep(0, scans)
    for (i in 1:numberofonsets) {
        onset <- onsets[i]
        for (j in onset:(onset + length - 1)) {
            stimulus[j] <- 1
        }
    }
    mygamma <- function(x, a1, a2, b1, b2, c) {
        d1 <- a1 * b1
        d2 <- a2 * b2
        c1 <- (x/d1)^a1
        c2 <- c * (x/d2)^a2
        res <- c1 * exp(-(x - d1)/b1) - c2 * exp(-(x - d2)/b2)
        res
    }
    hrf <- convolve(stimulus, mygamma(scans:1, a1, a2, b1/rt, 
        b2/rt, cc))
    dim(hrf) <- c(scans, 1)
    if (mean) {
        hrf - mean(hrf)
    }
    else {
        hrf
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
