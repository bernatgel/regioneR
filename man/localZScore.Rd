\name{localZScore}
\alias{localZScore}
\title{Local z-score}
\usage{
localZScore(A, pt, window, step, ...)
}
\arguments{
  \item{A}{a region set in any of the formats accepted by
  \code{\link{toGRanges}} (\code{\link{GenomicRanges}},
  \code{\link{data.frame}}, etc...)}

  \item{pt}{a permTestResult object}

  \item{window}{a window in wich the local Z-score will be
  calculated (bp)}

  \item{step}{the number of bp that divide each Z-score
  evaluation}

  \item{...}{further arguments to be passed to other
  methods.}
}
\value{
It returns a local z-score object
}
\description{
Evaluates tthe variation of the z-score in the vicinty of
the original region set
}
\examples{
genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
A <- createRandomRegions(nregions=20, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE)
B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))

pt <- overlapPermTest(A=A, B=B, ntimes=10, genome=genome, non.overlapping=FALSE)
plot(pt)

lz <- localZScore(A=A, B=B, pt=pt)
plot(lz)


pt2 <- permTest(A=A, B=B, ntimes=10, randomize.function=randomizeRegions, evaluate.function=list(overlap=numOverlaps, distance=meanDistance), genome=genome, non.overlapping=FALSE)
plot(pt2)

lz2 <- localZScore(A=A, B=B, pt2)
plot(lz2)
}
\seealso{
\code{\link{overlapPermTest}}, \code{\link{permTest}}
}

