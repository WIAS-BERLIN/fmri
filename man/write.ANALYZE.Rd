\name{write.ANALYZE}
\alias{write.ANALYZE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Write out 4 dimensional datacube in ANALYZE file format.
}
\usage{
write.ANALYZE(ttt, name = "data", size = "int", voxelsize = c(2, 2, 2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttt}{ 4 dimensional datacube }
  \item{name}{ name of the file to be written }
  \item{size}{ string. describing datatype (default "int") }
  \item{voxelsize}{ vector of length 3 describing the voxelsize in mm. }
}
\details{
  Writes datacube \code{ttt} to file named \code{name} in ANLYZE file
  format. \code{size} describes datatype and \code{voxelsize} the
  voxelsize in mm
}
\value{
  Nothing is returned.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{ This function requires the AnalyzeFMRI package and the
  f.write.analyze() function therein. This is because until now, we did
  not implement an output routine ourselves. }
\seealso{ \code{\link{read.ANALYZE}} }
\examples{
## Example 1
write.ANALYZE(as.integer(array(65526*runif(10*10*10*20),c(10,10,10,20))),"random_data")
}
\keyword{ SPM }% at least one, from doc/KEYWORDS
\keyword{ ANALYZE file format }% __ONLY ONE__ keyword per line
