\name{read.ANALYZE}
\alias{read.ANALYZE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.ANALYZE(prefix = "", numbered = FALSE, postfix = "", picstart = 1, numbpic = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prefix}{ ~~Describe \code{prefix} here~~ }
  \item{numbered}{ ~~Describe \code{numbered} here~~ }
  \item{postfix}{ ~~Describe \code{postfix} here~~ }
  \item{picstart}{ ~~Describe \code{picstart} here~~ }
  \item{numbpic}{ ~~Describe \code{numbpic} here~~ }
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
function (prefix = "", numbered = FALSE, postfix = "", picstart = 1, 
    numbpic = 1) 
{
    counter <- c(paste("00", 1:9, sep = ""), paste("0", 10:99, 
        sep = ""), paste(100:999, sep = ""))
    if (require(AnalyzeFMRI)) {
        if (numbered) {
            filename <- paste(prefix, counter[picstart], postfix, 
                ".img", sep = "")
        }
        else {
            filename <- paste(prefix, ".img", sep = "")
        }
        if (length(system(paste("ls", filename), TRUE, TRUE)) != 
            0) {
            ttt <- f.read.analyze.volume(filename)
            dt <- dim(ttt)
            cat(".")
            header <- f.read.analyze.header(filename)
            if (numbpic > 1) {
                for (i in (picstart + 1):(picstart + numbpic - 
                  1)) {
                  if (numbered) {
                    filename <- paste(prefix, counter[i], postfix, 
                      ".img", sep = "")
                  }
                  else {
                    filename <- paste(prefix, ".img", sep = "")
                  }
                  a <- f.read.analyze.volume(filename)
                  if (sum() != 0) 
                    cat("Error: wrong spatial dimension in picture", 
                      i)
                  ttt <- c(ttt, a)
                  dt[4] <- dt[4] + dim(a)[4]
                  cat(".")
                }
            }
            cat("\n")
            dim(ttt) <- c(dt)
            invisible(list(ttt = ttt, header = header))
        }
        else {
            warning(paste("Error: file", filename, "does not exist!"))
            list(ttt = NULL, header = NULL)
        }
    }
    else {
        warning("Error: library AnalyzeFMRI not found\n")
        list(ttt = NULL, header = NULL)
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
