\name{overlapPermTest}
\alias{overlapPermTest}
\title{Permutation Test for Overlap}
\usage{
overlapPermTest (A, B, alternative="auto", ...)
}
\arguments{
  \item{A}{a region set in any of the accepted formats by
  \code{\link{toGRanges}} (\code{\link{GenomicRanges}},
  \code{\link{data.frame}}, etc...)}

  \item{B}{a region set in any of the accepted formats by
  \code{\link{toGRanges}} (\code{\link{GenomicRanges}},
  \code{\link{data.frame}}, etc...)}

  \item{alternative}{the alternative hypothesis must be one
  of \code{"greater"}, \code{"less"} or \code{"auto"}. If
  \code{"auto"}, the alternative will be decided depending
  on the data.}

  \item{...}{further arguments to be passed to or from
  methods.}
}
\value{
A list of class \code{permTestResults} containing the
following components: \itemize{ \item \bold{\code{pval}}
the p-value of the test. \item \bold{\code{ntimes}} the
number of permutations. \item \bold{\code{alternative}} a
character string describing the alternative hypotesis.
\item \bold{\code{observed}} the value of the statistic for
the original data set. \item \bold{\code{permuted}} the
values of the statistic for each permuted data set. \item
\bold{\code{zscore}} the value of the standard score.
\code{(observed-\link{mean}(permuted))/\link{sd}(permuted)}
}
}
\description{
Performs a permutation test to see if there is an
association in overlap between a region set A and a region
set B creating random regions through the genome.
}
\examples{
genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE)
B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))

pt <- overlapPermTest(A=A, B=B, ntimes=10, genome=genome, non.overlapping=FALSE, verbose=TRUE)
summary(pt)
plot(pt)
plot(pt, plotType="Tailed")
}
\seealso{
\code{\link{overlapGraphicalSummary}},
\code{\link{overlapRegions}}, \code{\link{toDataframe}},
\code{\link{toGRanges}}, \code{\link{permTest}}
}

