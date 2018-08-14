#' Mean Distance
#' 
#' @description 
#' Computes the mean distance of regions in A to the nearest element in B
#' 
#' @usage 
#' meanDistance(A, B, ...)
#' 
#' @note 
#' If a region in A is in a chromosome where no B region is, it will be ignored and removed from the mean computation.
#' 
#' @param A   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param ... any additional parameter needed
#' 
#' @return
#' The mean of the distances of each region in A to the nearest region in B.
#' 
#' @examples
#' A <- data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40))
#' 
#' B <- data.frame("chr1", 25, 35)
#' 
#' meanDistance(A, B)
#' 
#' @export meanDistance
#' 
#' @importFrom GenomicRanges distanceToNearest


meanDistance <- function(A, B, ...) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")
  
  A <- toGRanges(A)
  B <- toGRanges(B)
  
  d <- GenomicRanges::distanceToNearest(A, B)
    
  return(mean(as.matrix(d@elementMetadata@listData$distance)[,1], na.rm=TRUE)) #--> BioC 2.13 
  #return(mean(d@listData$distance, na.rm=TRUE)) #--> BioC 2.11
  
}



