\name{numOverlaps}
\alias{numOverlaps}
\title{Number Of Overlaps}
\usage{
numOverlaps(A, B, count.once=FALSE, ...)
}
\arguments{
  \item{A}{a region set in any of the formats accepted by
  \code{\link{toGRanges}} (\code{\link{GenomicRanges}},
  \code{\link{data.frame}}, etc...)}

  \item{B}{a region set in any of the formats accepted by
  \code{\link{toGRanges}} (\code{\link{GenomicRanges}},
  \code{\link{data.frame}}, etc...)}

  \item{count.once}{boolean indicating whether the overlap
  of multiple B regions with a single A region should be
  counted once or multiple times}

  \item{...}{any additional parameters needed}
}
\value{
It returns a numeric value that is the number of regions in
A overlapping at least one region in B.
}
\description{
Returns the number of regions in A overlapping any region
in B
}
\examples{
genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE)
B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))

numOverlaps(A, B)
numOverlaps(A, B, count.once=TRUE)
}
\seealso{
\code{\link{overlapPermTest}}, \code{\link{permTest}}
}

