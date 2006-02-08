\name{read.AFNI}
\alias{read.AFNI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O function }
\description{
  Read BRIK file.
}
\usage{
read.AFNI(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ name of the file (without extension) }
}
\details{
  The function reads a BRIK/HEAD file.
}
\value{
  \item{ttt}{4 dimensional datacube}
  \item{header}{header information list}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{write.AFNI}} }
\examples{
\dontrun{afni <- read.AFNI("afnifile")}
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
