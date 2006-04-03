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
  Object of class \code{fmridata} with the following list entries:
  \item{ttt}{4 dimensional datacube}
  \item{header}{header information list}
  \item{format}{data source. string HEAD/BRIK}
  \item{delta}{voxel size in mm}
  \item{origin}{position of the datacube origin}
  \item{orient}{data orientation code. see AFNI documentation}
  \item{dim}{dimension of the datacube}
  \item{weights}{weights vector coding the different voxel sizes in x,
    y, z}
  \item{mask}{head mask}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\references{R. W. Cox (1996). AFNI: Software for analysis and visualization of
  functional magnetic resonance neuroimages. Computers and Biomed. Res. 
  29:162-173.
}
\seealso{ \code{\link{write.AFNI}} }
\examples{
\dontrun{afni <- read.AFNI("afnifile")}
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
