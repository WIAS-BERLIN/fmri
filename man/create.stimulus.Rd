\name{create.stimulus}
\alias{create.stimulus}

\title{Creates stimulus vector}
\description{Creates a vector of the convolution of the design function with
  hemodynamic response function given by Worsely et al. (2002)}

\usage{create.stimulus(scans = 1, onsets = c(1), length = 1, rt = 3, mean 
= TRUE)}

\arguments{
\item{scans}{Integer specifying the number of scans for the
  experiment. This will be the length of the resulting vector}
\item{onsets}{Vector of onset times of the stimulus.}
\item{length}{Length of the stimulus in scans starting at the onset times.}
\item{rt}{Repeat time between scans}
\item{mean}{Boolean variable indicating whether the resulting stimulus
  vector should have zero mean.}
}

\value{A vector of length \code{scans} containing the convolution of the
  stimulus with the hemodynamic response function.}

\examples{
hrf <- create.stimulus(107,c(22,52,82),11)
}
