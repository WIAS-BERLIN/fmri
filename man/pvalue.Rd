\name{pvalue}
\alias{pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pvalue(z, i, j, k, rx, ry, rz, type = "norm", df = 4, dim = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ ~~Describe \code{z} here~~ }
  \item{i}{ ~~Describe \code{i} here~~ }
  \item{j}{ ~~Describe \code{j} here~~ }
  \item{k}{ ~~Describe \code{k} here~~ }
  \item{rx}{ ~~Describe \code{rx} here~~ }
  \item{ry}{ ~~Describe \code{ry} here~~ }
  \item{rz}{ ~~Describe \code{rz} here~~ }
  \item{type}{ ~~Describe \code{type} here~~ }
  \item{df}{ ~~Describe \code{df} here~~ }
  \item{dim}{ ~~Describe \code{dim} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (z, i, j, k, rx, ry, rz, type = "norm", df = 4, dim = 3) 
{
    chihelp <- function(z, df, d) {
        z^((df - d)/2) * exp(-z/2)/2^((df - 2)/2)/gamma(df/2)
    }
    rho0 <- 0
    rho1 <- 0
    rho2 <- 0
    rho3 <- 0
    if (type == "norm") {
        rho0 <- 1 - pnorm(z)
        rho1 <- (4 * log(2))^0.5/(2 * pi) * exp(-z * z/2)
        rho2 <- (4 * log(2))/(2 * pi)^1.5 * exp(-z * z/2) * z
        rho3 <- (4 * log(2))^1.5/(2 * pi)^2 * exp(-z * z/2) * 
            (z * z - 1)
    }
    if (type == "t") {
        rho0 <- 1 - pt(z, df)
        rho1 <- (4 * log(2))^0.5/(2 * pi) * (1 + z * z/df)^(-0.5 * 
            (df - 1))
        rho2 <- (4 * log(2))/(2 * pi)^1.5 * gamma((df + 1)/2)/gamma(df/2)/(df/2)^0.5 * 
            (1 + z * z/df)^(-0.5 * (df - 1)) * z
        rho3 <- (4 * log(2))^1.5/(2 * pi)^2 * (1 + z * z/df)^(-0.5 * 
            (df - 1)) * ((df - 1)/df * z * z - 1)
    }
    if (type == "chisq") {
        rho0 <- 1 - pchisq(z, df)
        rho1 <- (4 * log(2))^0.5/(2 * pi)^0.5 * chihelp(z, df, 
            1)
        rho2 <- (4 * log(2))/(2 * pi) * chihelp(z, df, 2) * (z - 
            df + 1)
        rho3 <- (4 * log(2))^1.5/(2 * pi)^1.5 * chihelp(z, df, 
            3) * (z * z - (2 * df - 1) * z + (df - 1) * (df - 
            2))
    }
    r0 <- 1
    r1 <- (i - 1) * rx + (j - 1) * ry + (k - 1) * rz
    r2 <- (i - 1) * (j - 1) * rx * ry + (j - 1) * (k - 1) * ry * 
        rz + (i - 1) * (k - 1) * rx * rz
    r3 <- (i - 1) * (j - 1) * (k - 1) * rx * ry * rz
    if (dim == 2) {
        rho0 * r0 + rho1 * r1 + rho2 * r2
    }
    else {
        rho0 * r0 + rho1 * r1 + rho2 * r2 + rho3 * r3
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
