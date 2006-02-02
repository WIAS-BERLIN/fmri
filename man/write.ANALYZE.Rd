\name{write.ANALYZE}
\alias{write.ANALYZE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.ANALYZE(ttt, name = "data", size = "int", voxelsize = c(2, 2, 2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttt}{ ~~Describe \code{ttt} here~~ }
  \item{name}{ ~~Describe \code{name} here~~ }
  \item{size}{ ~~Describe \code{size} here~~ }
  \item{voxelsize}{ ~~Describe \code{voxelsize} here~~ }
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
function (ttt, name = "data", size = "int", voxelsize = c(2, 
    2, 2)) 
{
    if (require(AnalyzeFMRI)) {
        f.write.analyze(ttt, file = name, size = size, voxelsize)
    }
    else {
        warning("Error: library AnalyzeFMRI not found\n")
        invisible(NULL)
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
