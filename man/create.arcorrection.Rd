\name{create.arcorrection}
\alias{create.arcorrection}

\title{Creates an autocorrelation matrix}
\description{This function creates an autocorrelation matrix of a AR(1)
  type with paramter \code{rho}.
}

\usage{create.arcorrection(scans, rho = 0)}

\arguments{
  \item{scans}{Integer specifying the number of scans for the
  experiment.}
  \item{rho}{Autocorrelation parameter}
}

\value{The autocorreltation matrix of dimension \code{scans}.}

\examples{
ar <- create.arcorrection(107, 0.3)
}
