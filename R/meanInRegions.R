#' Mean In Regions
#' 
#' @description 
#' Returns the mean of a value defined by a region set over another set of regions.
#' 
#' @usage 
#' meanInRegions(A, x, col.name=NULL, ...)
#'  
#' @param A   a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param x   a region set in any of the accepted formats with an additional column with a value associated to every region. Regions in \code{x} can be points (single base regions).
#' @param col.name  character indicating the name of the column. If NULL and if a column with the name "value" exist, it will be used. The 4th column will be used otherwise (or the 5th if 4th is the strand). 
#' @param ... any additional parameter needed
#' 
#' @return
#' It returns a numeric value that is the weighted mean of "value" defined in \code{x} over the regions in \code{A}. That is, the mean of the value of all 
#' regions in \code{x} overlapping each region in \code{A} weighted according to the number of bases overlapping. 
#' 
#' @seealso \code{\link{permTest}} 
#' 
#' @examples
#' 
#'  A <- data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40))
#'  
#'  positions <- sample(1:40,30)
#'  
#'  x <- data.frame("chr1", positions, positions, rnorm(30,4,1))
#'  
#'  meanInRegions(A, x)
#'  
#'  x <- GRanges(seqnames=x[,1],ranges=IRanges(x[,2],end=x[,2]),mcols=x[,3])
#'  
#'  meanInRegions(A, x)
#'  
#'  @export meanInRegions
  



meanInRegions <- function(A, x, col.name=NULL, ...) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(x)) stop("x is missing")
  
  A <- toGRanges(A)
  x <- toGRanges(x)
  
  if(length(mcols(x))<1) {
    stop("x does not have a values column")
  }
  
  if(!is.null(col.name)) {
    value.col <-grep(col.name, names(mcols(x)))
    if(length(value.col) > 0) {
      value.col <- value.col[1]
    } else {
      value.col <- 1
    }
  } else {
    value.col <- 1
  }
  
  if(!is.numeric(mcols(x)[,value.col])) {
    stop("the values column in x is not numeric")
  } 

  value.col.name <- names(mcols(x))[value.col]
      
  over <- overlapRegions(A=A, B=x, colB=value.col, get.bases=TRUE) #TODO: If a strand column is present, col.B should be value.col + 4? REDO!

  if(length(over)==0) {
    warning("NA returned. There is no overlap between x and A.")
    return(NA)
  }
  
  total.value <- sum(as.numeric(over$ov.bases * over[,value.col.name]))  #Using as.numeric to escape possible integer overflows
  total.overlap <- sum(as.numeric(over$ov.bases))
  
  return(total.value/total.overlap)
  
}


