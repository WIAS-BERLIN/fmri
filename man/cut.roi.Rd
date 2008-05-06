\name{cut.roi}
\alias{cut.roi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cut.roi(data, xind = 1:data$dim[1], yind = 1:data$dim[2], zind = 1:data$dim[3], tind = 1:data$dim[4])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{xind}{ ~~Describe \code{xind} here~~ }
  \item{yind}{ ~~Describe \code{yind} here~~ }
  \item{zind}{ ~~Describe \code{zind} here~~ }
  \item{tind}{ ~~Describe \code{tind} here~~ }
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
function (data, xind = 1:data$dim[1], yind = 1:data$dim[2], zind = 1:data$dim[3], 
    tind = 1:data$dim[4]) 
{
    ttt <- extract.data(data)[xind, yind, zind, tind]
    data$ttt <- writeBin(as.numeric(ttt), raw(), 4)
    data$dim <- c(length(xind), length(yind), length(zind), length(tind))
    data$mask <- data$mask[xind, yind, zind]
    data$roixa <- (data$roixa:data$roixe)[xind[1]]
    data$roixe <- (data$roixa:data$roixe)[xind[length(xind)]]
    data$roiya <- (data$roiya:data$roiye)[yind[1]]
    data$roiye <- (data$roiya:data$roiye)[yind[length(yind)]]
    data$roiza <- (data$roiza:data$roize)[zind[1]]
    data$roize <- (data$roiza:data$roize)[zind[length(zind)]]
    data$roit <- data$roit[tind]
    invisible(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
