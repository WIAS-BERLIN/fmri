\name{write.AFNI}
\alias{write.AFNI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.AFNI(file, ttt, label, note = "", origin = c(0, 0, 0), delta = c(4, 4, 4), idcode = "WIAS_noid")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ ~~Describe \code{file} here~~ }
  \item{ttt}{ ~~Describe \code{ttt} here~~ }
  \item{label}{ ~~Describe \code{label} here~~ }
  \item{note}{ ~~Describe \code{note} here~~ }
  \item{origin}{ ~~Describe \code{origin} here~~ }
  \item{delta}{ ~~Describe \code{delta} here~~ }
  \item{idcode}{ ~~Describe \code{idcode} here~~ }
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
function (file, ttt, label, note = "", origin = c(0, 0, 0), delta = c(4, 
    4, 4), idcode = "WIAS_noid") 
{
    AFNIheaderpart <- function(type, name, value) {
        a <- "\n"
        a <- paste(a, "type = ", type, "\n", sep = "")
        a <- paste(a, "name = ", name, "\n", sep = "")
        if (regexpr("string", type) == 1) {
            value <- paste("'", value, "~", sep = "")
            a <- paste(a, "count = ", nchar(value) - 1, "\n", 
                sep = "")
            a <- paste(a, value, "\n", sep = "")
        }
        else {
            a <- paste(a, "count = ", length(value), "\n", sep = "")
            j <- 0
            while (j < length(value)) {
                left <- length(value) - j
                if (left > 4) 
                  left <- 5
                a <- paste(a, paste(value[(j + 1):(j + left)], 
                  collapse = "  "), "\n", sep = "  ")
                j <- j + 5
            }
        }
        a
    }
    conhead <- file(paste(file, ".HEAD", sep = ""), "w")
    writeChar(AFNIheaderpart("string-attribute", "HISTORY_NOTE", 
        note), conhead, eos = NULL)
    writeChar(AFNIheaderpart("string-attribute", "TYPESTRING", 
        "3DIM_HEAD_FUNC"), conhead, eos = NULL)
    writeChar(AFNIheaderpart("string-attribute", "IDCODE_STRING", 
        idcode), conhead, eos = NULL)
    writeChar(AFNIheaderpart("string-attribute", "IDCODE_DATE", 
        date()), conhead, eos = NULL)
    writeChar(AFNIheaderpart("integer-attribute", "SCENE_DATA", 
        c(0, 11, 1, -999, -999, -999, -999, -999)), conhead, 
        eos = NULL)
    writeChar(AFNIheaderpart("integer-attribute", "ORIENT_SPECIFIC", 
        c(0, 3, 4)), conhead, eos = NULL)
    writeChar(AFNIheaderpart("float-attribute", "ORIGIN", origin), 
        conhead, eos = NULL)
    writeChar(AFNIheaderpart("float-attribute", "DELTA", delta), 
        conhead, eos = NULL)
    minmax <- function(y) {
        r <- NULL
        for (k in 1:dim(y)[4]) {
            r <- c(r, min(y[, , , k]), max(y[, , , k]))
        }
        r
    }
    mm <- minmax(ttt)
    writeChar(AFNIheaderpart("float-attribute", "BRICK_STATS", 
        mm), conhead, eos = NULL)
    writeChar(AFNIheaderpart("integer-attribute", "DATASET_RANK", 
        c(3, dim(ttt)[4], 0, 0, 0, 0, 0, 0)), conhead, eos = NULL)
    writeChar(AFNIheaderpart("integer-attribute", "DATASET_DIMENSIONS", 
        c(dim(ttt)[1:3], 0, 0)), conhead, eos = NULL)
    writeChar(AFNIheaderpart("integer-attribute", "BRICK_TYPES", 
        rep(1, dim(ttt)[4])), conhead, eos = NULL)
    scale <- rep(0, dim(ttt)[4])
    for (k in 1:dim(ttt)[4]) {
        scale[k] <- max(abs(mm[2 * k - 1]), abs(mm[2 * k]))/32767
        ttt[, , , k] <- ttt[, , , k]/scale[k]
    }
    writeChar(AFNIheaderpart("float-attribute", "BRICK_FLOAT_FACS", 
        scale), conhead, eos = NULL)
    writeChar(AFNIheaderpart("string-attribute", "BRICK_LABS", 
        paste(label, collapse = "~")), conhead, eos = NULL)
    writeChar(AFNIheaderpart("string-attribute", "BYTEORDER_STRING", 
        "MSB_FIRST"), conhead, eos = NULL)
    close(conhead)
    conbrik <- file(paste(file, ".BRIK", sep = ""), "wb")
    dim(ttt) <- NULL
    writeBin(as.integer(ttt), conbrik, size = 2, endian = "big")
    close(conbrik)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
