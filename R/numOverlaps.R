#' Number Of Overlaps
#' 
#' @description 
#' Returns the number of regions in A overlapping any region in B
#' 
#' @usage 
#' numOverlaps(A, B, count.once=FALSE, ...)
#'  
#' @param A   a region set in any of the formats accepted by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B   a region set in any of the formats accepted by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param count.once boolean indicating whether the overlap of multiple B regions with a single A region should be counted once or multiple times
#' @param ... any additional parameters needed
#' 
#' @return
#' It returns a numeric value that is the number of regions in A overlapping at least one region in B.
#' 
#' @seealso \code{\link{overlapPermTest}}, \code{\link{permTest}}
#' 
#' @examples
#' 
#' genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
#' A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
#' B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
#' 
#' numOverlaps(A, B)
#' numOverlaps(A, B, count.once=TRUE)
#'  
#' @export numOverlaps
  


numOverlaps <- function(A, B, count.once=FALSE, ...) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")
  
  if(count.once) {
    return(length(which(overlapRegions(A, B, only.boolean=TRUE, ...))))
  } else {
    return(overlapRegions(A, B, only.count=TRUE, ...))  
  }
}


