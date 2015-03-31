#' Common Regions
#' 
#' @description 
#' Returns the regions that are common in two region sets, its intersection.
#' 
#' @note 
#' All metadata (additional columns in the region set in addition to chromosome, start and end) will be ignored and not present in the returned region set.
#' 
#' 
#' @usage 
#' commonRegions(A, B)
#' 
#' @param A   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' 
#' @return
#' It returns a \code{\link{GenomicRanges}} object with the regions present in both region sets.
#' 
#' @seealso \code{\link{plotRegions}}, \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{subtractRegions}}, \code{\link{splitRegions}}, \code{\link{extendRegions}}, \code{\link{joinRegions}}, \code{\link{mergeRegions}}, \code{\link{overlapRegions}}
#' 
#' @examples
#' A <- data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40))
#' 
#' B <- data.frame("chr1", 25, 35)
#' 
#' commons <- commonRegions(A, B)
#' 
#' plotRegions(list(A, B, commons), chromosome="chr1", regions.labels=c("A", "B", "common"), regions.colors=3:1)
#' 
#' @export commonRegions
#' 
#' 
#' 


commonRegions <- function(A, B) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")

  A <- toGRanges(A)
  B <- toGRanges(B)
  
  intersect(A, B)
  
}
