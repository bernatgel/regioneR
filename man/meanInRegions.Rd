\name{meanInRegions}
\alias{meanInRegions}
\title{Mean In Regions}
\usage{
meanInRegions(A, x, col.name=NULL, ...)
}
\arguments{
  \item{A}{a region set in any of the accepted formats by
  \code{\link{toGRanges}} (\code{\link{GenomicRanges}},
  \code{\link{data.frame}}, etc...)}

  \item{x}{a region set in any of the accepted formats with
  an additional column with a value associated to every
  region. Regions in \code{x} can be points (single base
  regions).}

  \item{col.name}{character indicating the name of the
  column. If NULL and if a column with the name "value"
  exist, it will be used. The 4th column will be used
  otherwise (or the 5th if 4th is the strand).}

  \item{...}{any additional parameter needed}
}
\value{
It returns a numeric value that is the weighted mean of
"value" defined in \code{x} over the regions in \code{A}.
That is, the mean of the value of all regions in \code{x}
overlapping each region in \code{A} weighted according to
the number of bases overlapping.
}
\description{
Returns the mean of a value defined by a region set over
another set of regions.
}
\examples{
A <- data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40))

 positions <- sample(1:40,30)

 x <- data.frame("chr1", positions, positions, rnorm(30,4,1))

 meanInRegions(A, x)

 x <- GRanges(seqnames=x[,1],ranges=IRanges(x[,2],end=x[,2]),mcols=x[,3])

 meanInRegions(A, x)
}
\seealso{
\code{\link{permTest}}
}

