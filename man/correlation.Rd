\name{correlation}
\alias{correlation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correlation(res, mask = array(1, dim = dim(res)[1:3]))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{ ~~Describe \code{res} here~~ }
  \item{mask}{ ~~Describe \code{mask} here~~ }
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
function (res, mask = array(1, dim = dim(res)[1:3])) 
{
    meanpos <- function(a) mean(a[a != 0])
    varpos <- function(a) var(a[a != 0])
    dr <- dim(res)
    if (length(dim(mask)) == 3) {
        if (sum(dim(mask)[1:3] != dr[1:3]) == 0) {
            mask <- rep(mask, dr[4])
            dim(mask) <- dr
            vrm <- varpos(res * mask)
            x <- meanpos(res[-1, , , ] * res[-dr[1], , , ] * 
                mask[-1, , , ])/vrm
            y <- meanpos(res[, -1, , ] * res[, -dr[2], , ] * 
                mask[, -1, , ])/vrm
            z <- meanpos(res[, , -1, ] * res[, , -dr[3], ] * 
                mask[, , -1, ])/vrm
            c(x, y, z)
        }
        else {
            warning("Error: dimension of mask and residui matrices do not match\n")
        }
    }
    else {
        warning("Error: can only handle 3 dimensional arrays\n")
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
