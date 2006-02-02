\name{gkernsm}
\alias{gkernsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gkernsm(y, h = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{h}{ ~~Describe \code{h} here~~ }
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
function (y, h = 1) 
{
    grid <- function(d) {
        d0 <- d\%/\%2 + 1
        gd <- seq(0, 1, length = d0)
        if (2 * d0 == d + 1) 
            gd <- c(gd, -gd[d0:2])
        else gd <- c(gd, -gd[(d0 - 1):2])
        gd
    }
    dy <- dim(y)
    if (is.null(dy)) 
        dy <- length(y)
    ldy <- length(dy)
    if (length(h) != ldy) 
        h <- rep(h[1], ldy)
    kern <- switch(ldy, dnorm(grid(dy), 0, 2 * h/dy), outer(dnorm(grid(dy[1]), 
        0, 2 * h[1]/dy[1]), dnorm(grid(dy[2]), 0, 2 * h[2]/dy[2]), 
        "*"), outer(outer(dnorm(grid(dy[1]), 0, 2 * h[1]/dy[1]), 
        dnorm(grid(dy[2]), 0, 2 * h[2]/dy[2]), "*"), dnorm(grid(dy[3]), 
        0, 2 * h[3]/dy[3]), "*"))
    kern <- kern/sum(kern)
    kernsq <- sum(kern^2)
    list(gkernsm = convolve(y, kern, conj = TRUE), kernsq = kernsq)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
