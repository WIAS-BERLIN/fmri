\name{read.ANALYZE}
\alias{read.ANALYZE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Read fMRI data from ANALYZE file.
}
\synopsis{
read.ANALYZE(prefix = "", numbered = FALSE, postfix = "", picstart = 1, numbpic = 1)
}
\usage{
read.ANALYZE(prefix = "run", numbered = TRUE, postfix = "_new", 
                                 picstart = 40, numbpic = 120)
read.ANALYZE(prefix = "filename")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prefix}{ string. part of the file name before the number }
  \item{numbered}{ logical. if \code{FALSE} only \code{prefix} is taken
    as file name.}
  \item{postfix}{ string. part of the file name after the number}
  \item{picstart}{ number of the first picture to be read. }
  \item{numbpic}{ number of pictures to be read }
}
\details{
  This function read in fMRI data file(s) in ANALYZE format as used in
  SPM. It takes the string in \code{prefix} and \code{postfix} and
  number of the form "007" in between to create the file name.

  If \code{numbered} is \code{FALSE}, only the string in \code{prefix}
  is used for file name (default).

  The number is assumed to be 3 digits (including leading zeros). First
  number is given is \code{picstart}, while \code{numbpic} defines the
  number of total numbur of pictures to be read. Data in multiple files
  will be combined into one four dimensional datacube.
}
\value{
  \item{ttt}{four dimensional data cube}
  \item{header}{header information of the data}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{ This function requires the AnalyzeFMRI package and the
  f.read.analyze.volume() and f.read.analyze.header() functions
  therein. This is because until now, we did not implement an input
  routine ourselves.

  Since numbering and naming of ANALYZE files widely varies, this
  function may not meet your personal needs. See Details for description.
}
\seealso{ \code{\link{write.ANALYZE}} }
\examples{ \dontrun{analyze <- read.ANALYZE("run",TRUE,"_new",40,120)}}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
