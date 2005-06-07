\name{calculate.lm}
\alias{calculate.lm}

\title{Calculate the linear model}
\description{This function performs the linear regression of the data to
  the design matrix.
}

\usage{calculate.lm(ttt,z,ar,vtype = "var")}

\arguments{
  \item{ttt}{Data matrix}
  \item{z}{Design matrix}
  \item{ar}{autocorrelation matrix}
  \item{vtype}{Type of determination for the varianz estimation}
}

\value{
  \item{beta}{Paramter estimate}
  \item{varianz}{Estimate of the variance of the parameters.}
  \item{residuals}{Residui}
}

\examples{
hrf <- create.stimulus(107, c(22, 52, 82), 12)
dim(hrf) <- c(107,1)
z <- create.designmatrix(hrf,2)
a <- create.arcorrection(107, rho=0.3)
ttt <- rnorm(1000)
dim(ttt) <- c(10,10,10)
model <- calculate.lm(ttt,z,a)
}
