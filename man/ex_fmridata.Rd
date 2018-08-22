\name{ex_fmridata}
\alias{ex_fmridata}
\docType{data}
\title{
Example fMRI data of the class of fmridata
}
\description{
Example fMRI data of the class of fmridata
}
\usage{data("ex_fmridata")}
\format{
  The format is:
List of 5
 $ ttt    : raw [1:47015800] c1 b3 ec c0 ...
 $ dim    : num [1:4] 65 65 26 107
 $ weights: num [1:3] 1 1 2
 $ mask   : num [1:65, 1:65, 1:26] 1 1 1 1 1 1 1 1 1 1 ...
 $ delta  : num [1:4] 1 1 1 1
 - attr(*, "class")= chr "fmridata"
}
\examples{
data(ex_fmridata)
## maybe str(ex_fmridata) ; plot(ex_fmridata) ...
}
\keyword{datasets}
