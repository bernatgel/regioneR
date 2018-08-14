#' Subtract Regions
#' 
#' @description 
#' Function for subtracting a region set from another region set.
#' 
#' @details
#' This function returns the regions in A minus the parts of them overlapping the regions in B. Overlapping regions in the result will be fused.
#' 
#' The implementation relies completely in the \code{setdiff} function from \code{IRanges} package.
#' 
#' @usage subtractRegions(A, B)
#' 
#' @param A   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' 
#' @return A GenomicRanges object
#' 
#' @examples 
#' A <- data.frame(chr=1, start=c(1, 15, 24, 31), end=c(10, 20, 30, 35))
#' 
#' B <- data.frame(chr=1, start=c(2, 12, 24, 35), end=c(5, 25, 29, 40))
#' 
#' subtract <- subtractRegions(A, B)
#' 
#' plotRegions(list(A, B, subtract), chromosome=1, regions.labels=c("A", "B", "subtract"), regions.colors=3:1)
#' 
#' @export subtractRegions
#' 


subtractRegions <- function(A, B) {

  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")
    
  A <- toGRanges(A)
  B <- toGRanges(B)
  
  if(length(A)==0 | length(B)==0) { return(A) }
  
  C <- GenomicRanges::setdiff(A, B) #Use the functionality available in GRanges
  
  return(C) 
    
}

