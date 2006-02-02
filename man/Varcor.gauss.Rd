\name{Varcor.gauss}
\alias{Varcor.gauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Varcor.gauss(h, interv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{ ~~Describe \code{h} here~~ }
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
function (h, interv = 1) 
{
    h <- h/2.3548 * interv
    ih <- trunc(4 * h) + 1
    dx <- 2 * ih + 1
    d <- length(h)
    penl <- dnorm(((-ih[1]):ih[1])/h[1])
    if (d == 2) 
        penl <- outer(penl, dnorm(((-ih[2]):ih[2])/h[2]), "*")
    if (d == 3) 
        penl <- outer(penl, outer(dnorm(((-ih[2]):ih[2])/h[2]), 
            dnorm(((-ih[3]):ih[3])/h[3]), "*"), "*")
    2 * sum(penl)^2/sum(diff(penl, interv)^2)/interv^(d)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
