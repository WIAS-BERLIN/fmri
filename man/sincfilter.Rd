\name{sincfilter}
\alias{sincfilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for sinc-interpolation
}
\description{
Performs sinc interpolation for a equidistant time series \code{x}
to times \code{t}.
}
\usage{
sincfilter(t, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
  vector of new time points
}
  \item{x}{
  observed time series at times \code{1:length(x)}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a vector of interpolated values of the time series at time points given in
\code{t}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
         Joerg Polzehl \email{polzehl@wias-berlin.de}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{slicetiming}}
}
\examples{
  x <- 1:107
  y <- rnorm(x)
  z <- sincfilter(seq(1,107,.01),y)
  plot(x, y, ylim=range(y,z))
  lines(seq(1,107,.01),z,col=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip  }% use one of  RShowDoc("KEYWORDS")
\keyword{utilities }% __ONLY ONE__ keyword per line
