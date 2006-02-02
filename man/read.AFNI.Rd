\name{read.AFNI}
\alias{read.AFNI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.AFNI(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ ~~Describe \code{file} here~~ }
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
function (file) 
{
    conhead <- file(paste(file, ".HEAD", sep = ""), "r")
    header <- readLines(conhead)
    close(conhead)
    types <- NULL
    args <- NULL
    counts <- NULL
    values <- NULL
    for (i in 1:length(header)) {
        if (regexpr("^type *= *", header[i]) != -1) {
            tmptype <- strsplit(header[i], " *= *")[[1]][2]
            types <- c(types, tmptype)
            args <- c(args, strsplit(header[i + 1], " *= *")[[1]][2])
            tmpcounts <- as.numeric(strsplit(header[i + 2], " *= *")[[1]][2])
            counts <- c(counts, tmpcounts)
            i <- i + 3
            tmpvalue <- ""
            while ((regexpr("^$", header[i]) == -1) && (i <= 
                length(header))) {
                tmpvalue <- paste(tmpvalue, header[i])
                i <- i + 1
            }
            tmpvalue <- sub("^ +", "", tmpvalue)
            if ((tmptype == "integer-attribute") || (tmptype == 
                "float-attribute")) {
                tmpvalue <- as.numeric(strsplit(tmpvalue, " +")[[1]])
            }
            values <- c(values, list(value = tmpvalue))
        }
    }
    names(values) <- args
    dx <- values$DATASET_DIMENSIONS[1]
    dy <- values$DATASET_DIMENSIONS[2]
    dz <- values$DATASET_DIMENSIONS[3]
    dt <- values$DATASET_RANK[2]
    scale <- values$BRICK_FLOAT_FACS
    size <- file.info(paste(file, ".BRIK", sep = ""))$size/(dx * 
        dy * dz * dt)
    if (regexpr("MSB", values$BYTEORDER_STRING[1]) != -1) {
        endian <- "big"
    }
    else {
        endian <- "little"
    }
    if (as.integer(size) == size) {
        conbrik <- file(paste(file, ".BRIK", sep = ""), "rb")
        myttt <- readBin(conbrik, "int", n = dx * dy * dz * dt * 
            size, size = size, signed = TRUE, endian = endian)
        close(conbrik)
        dim(myttt) <- c(dx, dy, dz, dt)
        for (k in 1:dt) {
            if (scale[k] != 0) {
                cat("scale", k, "with", scale[k], "\n")
                cat(range(myttt[, , , k]), "\n")
                myttt[, , , k] <- scale[k] * myttt[, , , k]
                cat(range(myttt[, , , k]), "\n")
            }
        }
        list(ttt = myttt, header = values)
    }
    else {
        warning("Error reading file: Could not detect size per voxel\n")
        list(ttt = NULL, header = values)
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
