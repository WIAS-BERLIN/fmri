\name{write.ANALYZE}
\alias{write.ANALYZE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Write out 4 dimensional datacube in ANALYZE file format.
}
\usage{
write.ANALYZE(ttt, header=NULL, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttt}{ 4 dimensional datacube }
  \item{header}{ header information }
  \item{file}{ file name }
}
\details{
  Writes datacube \code{ttt} to file named \code{file} in ANLYZE file
  format. \code{header} is a list that contains the header information
  as documented by the Mayo Foundation. We give here a short summary. If
  a value is not provided, it will be tried to fill it with reasonable
  defaults, but do not expect fine results, if the entry has a special
  important meaning (h.i. pixdim). "datatype1" - 10 byte character,
  "dbname" - 18 byte character, "extents" - integer, "sessionerror" -
  integer, "regular" -character, "hkey"- character, "dimension" - 8
  integers, starting with the number of dimensions!, "unused" - 7
  integers, "datatype" - integer describing datatype usually "4",
  "bitpix" - integer, "dimun0" - integer, "pixdim" - 8 float starting
  with voxel dimensions, "voxoffset" - float, "funused" - 3 floats,
  "calmax" - float, "calmin" - float, "compressed" - float, "verified" -
  float, "glmax" - integer, "glmin" - integer, "describ" - 80 byte
  character, "auxfile"- 24 byte character, "orient" - character,
  "originator" - 5 integers (cf. ANALZYE docu!), "generated" - 10 byte
  character, "scannum" - 10 byte character, "patientid" - 10 byte
  character, "expdate" - 10 byte character, "exptime" - 10 byte
  character, "histun0" - 3 byte character, "views" - integer, "voladded"
  - integer, "startfield" - integer, "fieldskip" - integer, "omax" -
  integer, "omin" - integer, "smax" - integer, "smin" - integer       

}
\value{
  Nothing is returned.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{read.ANALYZE}}, \code{\link{write.AFNI}} }
\examples{
## Example 1
write.ANALYZE(array(as.integer(65526*runif(10*10*10*20)),c(10,10,10,20)),
              file="analyzefile")
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
