\name{Spatialvar.gauss}
\alias{Spatialvar.gauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Spatialvar.gauss(h, h0, d, interv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{ ~~Describe \code{h} here~~ }
  \item{h0}{ ~~Describe \code{h0} here~~ }
  \item{d}{ ~~Describe \code{d} here~~ }
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
function (h, h0, d, interv = 1) 
{
    h0 <- pmax(h0, 1e-05)
    h <- pmax(h, 1e-05)
    h <- h/2.3548 * interv
    if (length(h) == 1) 
        h <- rep(h, d)
    ih <- trunc(4 * h)
    ih <- pmax(1, ih)
    dx <- 2 * ih + 1
    penl <- dnorm(((-ih[1]):ih[1])/h[1])
    if (d == 2) 
        penl <- outer(dnorm(((-ih[1]):ih[1])/h[1]), dnorm(((-ih[2]):ih[2])/h[2]), 
            "*")
    if (d == 3) 
        penl <- outer(dnorm(((-ih[1]):ih[1])/h[1]), outer(dnorm(((-ih[2]):ih[2])/h[2]), 
            dnorm(((-ih[3]):ih[3])/h[3]), "*"), "*")
    dim(penl) <- dx
    h0 <- h0/2.3548 * interv
    if (length(h0) == 1) 
        h0 <- rep(h0, d)
    ih <- trunc(4 * h0)
    ih <- pmax(1, ih)
    dx0 <- 2 * ih + 1
    x <- ((-ih[1]):ih[1])/h0[1]
    penl0 <- dnorm(((-ih[1]):ih[1])/h0[1])
    if (d == 2) 
        penl0 <- outer(dnorm(((-ih[1]):ih[1])/h0[1]), dnorm(((-ih[2]):ih[2])/h0[2]), 
            "*")
    if (d == 3) 
        penl0 <- outer(dnorm(((-ih[1]):ih[1])/h0[1]), outer(dnorm(((-ih[2]):ih[2])/h0[2]), 
            dnorm(((-ih[3]):ih[3])/h0[3]), "*"), "*")
    dim(penl0) <- dx0
    penl0 <- penl0/sum(penl0)
    dz <- dx + dx0 - 1
    z <- array(0, dz)
    if (d == 1) {
        for (i1 in 1:dx0) {
            ind1 <- c(0:(i1 - 1), (dz - dx0 + i1):dz + 1)
            ind1 <- ind1[ind1 <= dz][-1]
            z[-ind1] <- z[-ind1] + penl * penl0[i1]
        }
    }
    else if (d == 2) {
        for (i1 in 1:dx0[1]) for (i2 in 1:dx0[2]) {
            ind1 <- c(0:(i1 - 1), (dz[1] - dx0[1] + i1):dz[1] + 
                1)
            ind1 <- ind1[ind1 <= dz[1]][-1]
            ind2 <- c(0:(i2 - 1), (dz[2] - dx0[2] + i2):dz[2] + 
                1)
            ind2 <- ind2[ind2 <= dz[2]][-1]
            z[-ind1, -ind2] <- z[-ind1, -ind2] + penl * penl0[i1, 
                i2]
        }
    }
    else if (d == 3) {
        for (i1 in 1:dx0[1]) for (i2 in 1:dx0[2]) for (i3 in 1:dx0[3]) {
            ind1 <- c(0:(i1 - 1), (dz[1] - dx0[1] + i1):dz[1] + 
                1)
            ind1 <- ind1[ind1 <= dz[1]][-1]
            ind2 <- c(0:(i2 - 1), (dz[2] - dx0[2] + i2):dz[2] + 
                1)
            ind2 <- ind2[ind2 <= dz[2]][-1]
            ind3 <- c(0:(i3 - 1), (dz[3] - dx0[3] + i3):dz[3] + 
                1)
            ind3 <- ind3[ind3 <= dz[3]][-1]
            z[-ind1, -ind2, -ind3] <- z[-ind1, -ind2, -ind3] + 
                penl * penl0[i1, i2, i3]
        }
    }
    sum(z^2)/sum(z)^2 * interv^d
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
