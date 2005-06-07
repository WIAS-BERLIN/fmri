\name{plot.fmri}
\alias{plot.fmri}

\title{Plot the activation areas}
\description{The function plots the result of the fmri analysis to
  screen or file. This is done slice by slice.
}

\usage{plot.fmri(signal, mask, anatomic, zlim=0, device="X11", file="plot.png")}

\arguments{
  \item{signal}{Signal strength}
  \item{mask}{Mask where signal has been found
  \item{anatomic}{Anatomic picture of same dimension}
  \item{zlim}{}
  \item{device}{Output device. Choose "png" for file output}
  \item{file}{Name of the file}
}

\value{Number of the graphic device opened. Close the device for writing
to disk in case of output to file.}

\examples{
ttt <- rnorm(1000)
dim(ttt) <- c(10,10,10)
mask <- array(1,dim=dim(ttt))
mask[ttt<0.5] <- 0
graph <- plot.fmri(mask,mask,ttt,device="png",file="mask.png")
dev.off(graph)
}
