#' Permutation Test for Overlap
#' 
#' @description
#' Performs a permutation test to see if the overlap between two sets of regions
#' A and B is higher (or lower) than expected by chance. It will internally
#' call \code{\link{permTest}} with the appropiate parameters to perform the 
#' permutation test. If B is a list or a GRangesList, it will perform one
#' permutation test per element of the list, testing the overlap between
#' A and each element of B independently.
#' 
#' @note \bold{IMPORTANT:} Since it uses \code{link{permTest}} internally, it
#' is possible to use most of the parameters of that function in 
#' \code{overlapPermTest}, including: \code{ntimes}, \code{force.parallel},
#' \code{min.parallel} and \code{verbose}. In addition, this function
#' accepts most parameters of the \code{\link{randomizeRegions}} function 
#' including \code{genome}, \code{mask}, \code{allow.overlaps} and 
#' \code{per.chromosome} and the parameters of \code{\link{numOverlaps}} such
#' as \code{count.once}.
#' 
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
#'  
#' C <- c(B, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
#' pt <- overlapPermTest(A=A, B=list(B=B, C=C), ntimes=10, genome=genome, non.overlapping=FALSE, verbose=TRUE)
#' summary(pt)
#' plot(pt)
#'  
#' @export overlapPermTest

#Convenience function to perform a a permutation test to assess the overlap between two different sets of regions A and B

overlapPermTest <- function(A, B, alternative="auto", ...) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")
  alternative <- match.arg(alternative,c("less","greater", "auto"))
  
  A <- toGRanges(A)
  
  if(methods::is(B, "GRangesList")) {
    B <- as.list(B)
  }
  
  #If there are multiple B's, create a list of curried functions 
  #(with parameter B pre-applied) and use that as the evaluation functions
  if(is.list(B) || is.vector(B)) {
    func.names <- NULL
    if(is.null(names(B))) {
      #if it's a vector of characters use them
      if(is.vector(B) && all(is.character(B))) {
        func.names <- B
      }
      #if it's a list of characters, use them
      if(is.list(B) && all(unlist(lapply(B, is.character)))) {
        func.names <- unlist(B)
      }
    }
    functs <- createFunctionsList(numOverlaps, param.name = "B", values = B, func.names = func.names)
    
    return(permTest(A=A, randomize.function=randomizeRegions, evaluate.function=functs, alternative=alternative, ...))
  } else {
    B <- toGRanges(B)  
    return(permTest(A=A, B=B, randomize.function=randomizeRegions, evaluate.function=numOverlaps, alternative=alternative, ...))    
  }
}
