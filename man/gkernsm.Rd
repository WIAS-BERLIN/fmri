\name{gkernsm}
\alias{gkernsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothness }
\description{
  Simple gaussian smoothing.
}
\usage{
gkernsm(y, h = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Data (1, 2, 3 dimensional)}
  \item{h}{ Bandwidth }
}
\details{
  Smooth the data with a gaussian kernel with andwidth \code{h}.
}
\value{
  \item{gkernsm}{Smoothed array.}
  \item{kernsq}{sum of squared kernel}
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\examples{
  # Example 1
  y <- array(rnorm(100*100),c(100,100))
  ysmooth <- gkernsm(y,5)$gkernsm
  oldpar <- par(mfrow=c(1,2))
  image(y)
  image(ysmooth)
  par(oldpar)
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
