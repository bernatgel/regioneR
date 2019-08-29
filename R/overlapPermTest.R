#' Permutation Test for Overlap
#' 
#' @description
#' Performs a permutation test to see if there is an association in overlap between a region set A and a region set B creating random regions through the genome.
#' 
#' @usage overlapPermTest (A, B, alternative="auto", ...)
#' 
#' @param A a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param alternative the alternative hypothesis must be one of \code{"greater"}, \code{"less"} or \code{"auto"}. If \code{"auto"}, the alternative will be decided depending on the data.
#' @param ... further arguments to be passed to or from methods.
#' 
#' @return
#' A list of class \code{permTestResults} containing the following components:
#' \itemize{
#' \item \bold{\code{pval}} the p-value of the test.
#' \item \bold{\code{ntimes}} the number of permutations.
#' \item \bold{\code{alternative}} a character string describing the alternative hypotesis.
#' \item \bold{\code{observed}} the value of the statistic for the original data set.
#' \item \bold{\code{permuted}} the values of the statistic for each permuted data set.
#' \item \bold{\code{zscore}} the value of the standard score. \code{(observed-\link{mean}(permuted))/\link{sd}(permuted)}
#' }
#' 
#' @seealso \code{\link{overlapGraphicalSummary}}, \code{\link{overlapRegions}}, \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{permTest}}
#' 
#' @examples
#' genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
#' A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
#' B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
#' 
#' pt <- overlapPermTest(A=A, B=B, ntimes=10, genome=genome, non.overlapping=FALSE, verbose=TRUE)
#' summary(pt)
#' plot(pt)
#' plot(pt, plotType="Tailed")  
#'  
#' @export overlapPermTest

#Convenience function to perform a a permutation test to assess the relation between two different sets of regions: A and B

overlapPermTest <- function(A, B, alternative="auto", ...) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")
  alternative<-match.arg(alternative,c("less","greater", "auto"))
  
  
  B <- toGRanges(B)
  A <- toGRanges(A)
  
  return(permTest(A=A, B=B, randomize.function=randomizeRegions, evaluate.function=numOverlaps, alternative=alternative, ...))
  
}
