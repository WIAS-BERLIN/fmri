\name{plot.signal}
\alias{plot.signal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.signal(signal, anatomic = array(0, dim = dim(signal)), pos = c(-1, -1, -1), device = "X11", file = "plot.png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal}{ ~~Describe \code{signal} here~~ }
  \item{anatomic}{ ~~Describe \code{anatomic} here~~ }
  \item{pos}{ ~~Describe \code{pos} here~~ }
  \item{device}{ ~~Describe \code{device} here~~ }
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
function (signal, anatomic = array(0, dim = dim(signal)), pos = c(-1, 
    -1, -1), device = "X11", file = "plot.png") 
{
    partition <- ceiling(sqrt(dim(anatomic)[3]))
    switch(device, png = png(filename = file, width = 200 * partition, 
        height = 200 * partition, pointsize = 12, bg = "transparent", 
        res = NA), jpeg = jpeg(filename = file, width = 200 * 
        partition, height = 200 * partition, quality = 100, pointsize = 12, 
        bg = "transparent", res = NA), ppm = bitmap(file, type = "ppm", 
        height = 2 * partition, width = 2 * partition, res = 64, 
        pointsize = 12), X11())
    oldpar <- par(mfrow = c(partition, partition), mar = c(0.25, 
        0.25, 0.25, 0.25), mgp = c(2, 1, 0))
    alim <- range(anatomic)
    zlim <- range(signal)
    for (i in 1:dim(anatomic)[3]) {
        image(anatomic[, , i], xaxt = "n", yaxt = "n", zlim = alim, 
            col = grey(0:255/255))
        if (any(signal[, , i])) 
            image(signal[, , i], zlim = zlim, col = c(0, heat.colors(512)), 
                add = TRUE)
        if (i == pos[3]) {
            lines(c(0, 1), c(pos[2], pos[2])/dim(anatomic)[2], 
                col = 2)
            lines(c(pos[1], pos[2])/dim(anatomic)[1], c(0, 1), 
                col = 2)
        }
    }
    par(oldpar)
    switch(device, png = dev.off(), jpeg = dev.off(), ppm = dev.off())
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
