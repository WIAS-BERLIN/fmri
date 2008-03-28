\name{fmriica}
\alias{fmriica}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Independent Component Analysis }
\description{
  The function performs Independent Component Analysis.
}
\usage{
fmriica(data, m = 3, method = "spatial", xind = NULL, yind = NULL, zind = NULL, tind = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{  Observation matrix (dimension Nxd) }
  \item{m}{ Number of independent components. }
  \item{method}{ ~~Describe \code{method} here~~ }
  \item{xind}{ ~~Describe \code{xind} here~~ }
  \item{yind}{ ~~Describe \code{yind} here~~ }
  \item{zind}{ ~~Describe \code{zind} here~~ }
  \item{tind}{ ~~Describe \code{tind} here~~ }
  \item{\dots}{ further arguments to fastICA }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{The function returns a list with components
    \item{ihat }{ Matrix containing the first m NGCA directions as columns.}
    \item{sdev }{ Standard deviations of the principal components of 
                 the thresholded ICA directions}
    \item{xhat }{first m components of the rotated data}
    \item{v }{If \code{keepv==TRUE} the set of directions \code{v^{(k)}}}
    \item{normv }{If \code{keepv==TRUE} the norm of each \code{v^{(k)}}.}
}
\references{ ~put references to the literature/web site here ~ }
\author{ J\"org Polzehl polzehl@wias-berlin.de  }
\seealso{ \code{\link{ngca}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }
