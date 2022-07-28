#' Resample Regions
#' 
#' @description 
#' Function for sampling a region set from a universe of region sets.
#' 
#' @usage 
#' resampleRegions(A, universe, per.chromosome=FALSE, ...)
#' 
#' @param A a region set in any of the formats accepted by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param universe a region set in any of the formats accepted by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param per.chromosome boolean indicating if sample must be by chromosome. 
#' @param ... further arguments to be passed to or from methods.
#' 
#' @return  a \code{\link{GenomicRanges}} object. A sample from the \code{univers} with the same length as A.
#' 
#' @seealso  \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{randomizeRegions}}, \code{\link{createRandomRegions}}
#' 
#' @examples 
#' universe <- data.frame(chr=1, start=c(1,15,24,40,50), end=c(10,20,30,45,55))
#' 
#' A <- data.frame(chr=1, start=c(2,12,28,35), end=c(5,25,33,43))
#' 
#' resampleRegions(A, universe, per.chromosome=TRUE)
#'  
#' @export resampleRegions
#' 



resampleRegions <- function(A, universe, per.chromosome=FALSE, ...) { 
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(universe)) stop("universe is missing")
  if(!is.logical(per.chromosome)) stop("per.chromosome must be logical")
  
  
  A <- toGRanges(A)
  universe <- toGRanges(universe)
  
  
  if(per.chromosome){
    chrResample <- function(chr) {
      Achr <- A[seqnames(A) == chr]
      universe.chr <- universe[seqnames(universe) == chr]
      resample.chr <- universe.chr[sample(1:length(universe.chr), length(Achr))]
      return(resample.chr)
    }
    
    chr.resampled <- lapply(as.list(seqlevels(A)), chrResample)
    resampled <- do.call(c, chr.resampled)
    
  }else{
    resampled <- universe[sample(1:length(universe), length(A))]
  }
  
  return(resampled)
  
}





