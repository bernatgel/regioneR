#' Extend Regions
#' 
#' @description 
#' Extends the regions a number of bases at each end. Negative numbers will reduce the region instead of enlarging it.
#' 
#' @note
#' If negative values are provided and the new extremes are "flipped", the function will fail. It does not check if the extended regions fit into the genome.
#' 
#' @usage 
#' extendRegions(A, extend.start=0, extend.end=0)
#' 
#' @param A a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param extend.start an integer. The number of bases to be subtracted from the start of the region.
#' @param extend.end an integer. The number of bases to be added at the end of the region.
#' 
#' @return
#' a \code{\link{GenomicRanges}} object with the extended regions.
#' 
#' @seealso \code{\link{plotRegions}}, \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{subtractRegions}}, \code{\link{splitRegions}}, \code{\link{overlapRegions}}, \code{\link{commonRegions}}, \code{\link{mergeRegions}}, \code{\link{joinRegions}}
#' 
#' @examples
#' A <- data.frame("chr1", c(10, 20, 30), c(13, 28, 40))
#' 
#' extend1 <- extendRegions(A, extend.start=5, extend.end=2)
#' 
#' extend2 <- extendRegions(A, extend.start=15)
#' 
#' extend3 <- extendRegions(A, extend.start=-1)
#' 
#' plotRegions(list(A, extend1, extend2, extend3), chromosome="chr1", regions.labels=c("A", "extend1", "extend2", "extend3"), regions.colors=4:1)
#' 
#' 
#' @export extendRegions


#It does'nt check for inclusion in the genome. Start may end up being less than 0 and end greater than chromosome length.
extendRegions <- function(A, extend.start=0, extend.end=0) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!is.numeric(extend.start)) stop("extend.start must be numeric")
  if(!is.numeric(extend.end)) stop("extend.end must be numeric")

  A <- toGRanges(A)
  
  start(A) <- start(A) - extend.start
  end(A) <- end(A) + extend.end
  
  return(A)
  
}