#' Unique Regions
#' 
#' @description 
#' Returns the regions unique to only one of the two region sets, that is, all parts of the genome covered by only one of the two region sets.
#' 
#' @note 
#' All metadata (additional columns in the region set in addition to chromosome, start and end) will be ignored and not present in the returned region set.
#' 
#' 
#' @usage 
#' uniqueRegions(A, B)
#' 
#' @param A   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' 
#' @return
#' It returns a \code{\link{GenomicRanges}} object with the regions unique to one of the region sets.
#' 
#' @seealso \code{\link{toGRanges}}, \code{\link{subtractRegions}}, \code{\link{commonRegions}}, \code{\link{mergeRegions}}
#' 
#' @examples
#' A <- data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40))
#' 
#' B <- data.frame("chr1", 25, 35)
#' 
#' uniques <- uniqueRegions(A, B)
#' 
#' plotRegions(list(A, B, uniques), chromosome="chr1", regions.labels=c("A", "B", "uniques"), regions.colors=3:1)
#' 
#' @export uniqueRegions
#' 


#The implementation is based on diff=union - intersection
uniqueRegions <- function(A, B) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")

  A <- toGRanges(A)
  B <- toGRanges(B)
  
  merged <- mergeRegions(A, B)
  common <- commonRegions(A, B)
  
  return(subtractRegions(merged,common))
  
}