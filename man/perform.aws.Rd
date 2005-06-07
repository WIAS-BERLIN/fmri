\name{perform.aws}
\alias{perform.aws}

\title{Perform adaptive weights smoothing}
\description{This functions performs adaptive weights smoothing on the
  parameter \code{beta}.
}

\usage{perform.aws(beta, varianz, hmax = 4, hinit = 1, weights =
c(1,1,1), qlambda = NULL)}

\arguments{
  \item{beta}{Data to smooth}
  \item{varianz}{Variance estimate for heteroscedastic smoothing}
  \item{hmax}{Paramter \code{hmax} of AWS procedure}
  \item{hinit}{Paramter \code{hinit} of AWS procedure}
  \item{weights}{Paramter \code{weights} of AWS procedure}
  \item{qlambda}{Paramter \code{qlambda} of AWS procedure}
}

\value{
  \item{hat}{Estimate of parameter after smoothing}
  \item{var}{Variance of the estimate}
}

\examples{
hrf <- create.stimulus(107, c(22, 52, 82), 12)
dim(hrf) <- c(107,1)
z <- create.designmatrix(hrf,2)
a <- create.arcorrection(107, rho=0.3)
ttt <- rnorm(1000)
dim(ttt) <- c(10,10,10)
model <- calculate.lm(ttt,z,a)
result4 <- analyze.fmri(model$beta[,,,1],model$varianz,hmax=4,weights=c(1,1,1.4))
}
