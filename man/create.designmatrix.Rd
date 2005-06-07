\name{create.designmatrix}
\alias{create.designmatrix}

\title{Creates a design matrix}
\description{This function creates a design matrix containing the
  stimuli given in \code{hrf} and a polynomial drift of order \code{order}
}

\usage{create.designmatrix(hrf, order = 0)}

\arguments{
  \item{hrf}{Matrix containing the different stimuli as rows.}
  \item{order}{Order of the polynomial drift for regression.}
}

\value{A matrix of dimension \code{dim(hrf) + order + 1} containing the
  design matrix.}
\seealso{
  \code{\link{create.stimulus}}
}

\examples{
hrf <- create.stimulus(107,c(22,52,82),11)
z <- create.designmatrix(hrf,4)
}
