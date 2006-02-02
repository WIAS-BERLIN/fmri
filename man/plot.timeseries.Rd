\name{plot.timeseries}
\alias{plot.timeseries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.timeseries(fmri, fit, sigvec, tvec, thr, origin, device = "X11", file = "plot.png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fmri}{ ~~Describe \code{fmri} here~~ }
  \item{fit}{ ~~Describe \code{fit} here~~ }
  \item{sigvec}{ ~~Describe \code{sigvec} here~~ }
  \item{tvec}{ ~~Describe \code{tvec} here~~ }
  \item{thr}{ ~~Describe \code{thr} here~~ }
  \item{origin}{ ~~Describe \code{origin} here~~ }
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
function (fmri, fit, sigvec, tvec, thr, origin, device = "X11", 
    file = "plot.png") 
{
    partition <- dim(fmri)[2:1]
    switch(device, png = png(filename = file, width = 200 * partition[2], 
        height = 200 * partition[1], pointsize = 12, bg = "transparent", 
        res = NA), jpeg = jpeg(filename = file, width = 200 * 
        partition[2], height = 200 * partition[1], quality = 100, 
        pointsize = 12, bg = "transparent", res = NA), ppm = bitmap(file, 
        type = "ppm", height = 2 * partition[1], width = 2 * 
            partition[2], res = 64, pointsize = 12), X11())
    oldpar <- par(mfrow = partition, mar = c(0.25, 0.25, 0.25, 
        0.25), mgp = c(2, 1, 0))
    ylim <- range(fmri)
    for (j in 1:dim(fmri)[2]) {
        for (i in 1:dim(fmri)[1]) {
            if (tvec[i, j, 1] > thr) {
                plot(fmri[i, j, ], cex = 0.5, xaxt = "n", yaxt = "n", 
                  col = "red")
            }
            else {
                plot(fmri[i, j, ], cex = 0.5, xaxt = "n", yaxt = "n")
            }
            lines(fit[i, j, ], col = 3)
            text(10, max(fmri[i, j, ]), pos = 4, paste(origin[1] + 
                i - 1, origin[2] + j - 1, sep = "/"))
            text(10, min(fmri[i, j, ]), pos = 4, paste(paste(signif(sigvec[i, 
                j, ], digits = 2), signif(tvec[i, j, ], digits = 2), 
                sep = "/"), collapse = "---"))
        }
    }
    par(oldpar)
    switch(device, png = dev.off(), jpeg = dev.off(), ppm = dev.off())
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
