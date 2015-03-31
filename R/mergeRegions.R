#' Merge Regions
#' 
#' @description 
#' Merges the overlapping regions from two region sets. The two region sets are first merged into one and then overlapping regions are fused. 
#' 
#' @note 
#' All metadata (additional columns in the region set in addition to chromosome, start and end) will be ignored and not present in the returned region set.
#'
#' The implementation relies completely in the \code{\link{reduce}} function from \code{IRanges} package.
#'  
#' @usage 
#' mergeRegions(A, B)
#' 
#' @param A   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' 
#' @return
#' It returns a \code{\link{GenomicRanges}} object with the regions resulting from the merging process. Any two overlapping regions from any of the two sets will be fused into one.
#' 
#' @seealso \code{\link{plotRegions}}, \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{subtractRegions}}, \code{\link{splitRegions}}, \code{\link{extendRegions}}, \code{\link{joinRegions}}, \code{\link{commonRegions}}, \code{\link{overlapRegions}}
#' 
#' @examples
#' A <- data.frame("chr1", c(1, 5, 20, 30), c(8, 13, 28, 40), x=c(1,2,3,4), y=c("a", "b", "c", "d"))
#' 
#' B <- data.frame("chr1", 25, 35)
#' 
#' merges <- mergeRegions(A, B)
#' 
#' plotRegions(list(A, B, merges), chromosome="chr1", regions.labels=c("A", "B", "merges"), regions.colors=3:1)
#' 
#' @export mergeRegions


mergeRegions <- function(A, B) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")

  A <- toGRanges(A)
  B <- toGRanges(B)
  
  C <- c(A,B, ignore.mcols=TRUE)  
  
  return(reduce(C))
  
}