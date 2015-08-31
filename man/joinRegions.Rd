\name{joinRegions}
\alias{joinRegions}
\title{Join Regions}
\usage{
joinRegions(A, min.dist=1)
}
\arguments{
  \item{A}{a region set in any of the accepted formats by
  \code{\link{toGRanges}} (\code{\link{GenomicRanges}},
  \code{\link{data.frame}}, etc...)}

  \item{min.dist}{an integer indicating the minimum
  distance required between two regions in order to not
  fuse them. Any pair of regions closer than
  \code{min.dist} bases will be fused in a larger region.
  Defaults to 1, so it will only join overlapping regions.}
}
\value{
It returns a \code{\link{GenomicRanges}} object with the
regions resulting from the joining process.
}
\description{
Joins the regions from a region set A that are less than
\code{min.dist} bases apart.
}
\note{
All metadata (additional columns in the region set in
addition to chromosome, start and end) will be ignored and
not present in the returned region set.

The implementation relies completely in the
\code{\link{reduce}} function from \code{IRanges} package.
}
\examples{
A <- data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40))

join1 <- joinRegions(A)

join2 <- joinRegions(A, min.dist=3)

join3 <- joinRegions(A, min.dist=10)

plotRegions(list(A, join1, join2, join3), chromosome="chr1", regions.labels=c("A", "join1", "join2", "join3"), regions.colors=4:1)
}
\seealso{
\code{\link{plotRegions}}, \code{\link{toDataframe}},
\code{\link{toGRanges}}, \code{\link{subtractRegions}},
\code{\link{splitRegions}}, \code{\link{extendRegions}},
\code{\link{commonRegions}}, \code{\link{mergeRegions}},
\code{\link{overlapRegions}}
}

