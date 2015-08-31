\name{plotRegions}
\alias{plotRegions}
\title{Plot Regions}
\usage{
plotRegions(x, chromosome, start=NULL, end=NULL, regions.labels=NULL, regions.colors=NULL, ...)
}
\arguments{
  \item{x}{list of objects to be ploted.}

  \item{chromosome}{character or numeric value indicating
  which chromosome you want to plot.}

  \item{start}{numeric value indicating from which position
  you want to plot.}

  \item{end}{numeric value indicating to which position you
  want to plot.}

  \item{regions.labels}{vector indicating the labels for
  the y axes. It must have the same length as x.}

  \item{regions.colors}{character vector indicating the
  colors for the plotted regions. It must have the same
  length as x.}

  \item{...}{Arguments to be passed to methods, such as
  graphical parameters (see \code{\link{par}}).}
}
\value{
A plot is created on the current graphics device.
}
\description{
Plots sets of regions
}
\examples{
A <- data.frame(chr=1, start=c(1,15,24,40,50), end=c(10,20,30,45,55))

B <- data.frame(chr=1, start=c(2,12,28,35), end=c(5,25,33,43))

plotRegions(list(A,B), chromosome=1, regions.labels=c("A","B"), regions.colors=3:2)
}

