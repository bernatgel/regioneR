#' toDataframe
#' 
#' @description
#' Transforms a \code{\link{GRanges}} object or a \code{\link{data.frame}}containing a region set into a \code{\link{data.frame}}.
#' 
#' @details
#' If the oject is of class \code{\link{data.frame}}, it will be returned untouched.
#' 
#' @usage toDataframe(A, stranded=FALSE)
#' 
#' @param A a \code{\link{GRanges}} object.
#' @param stranded    (only used when A is a \code{\link{GRanges}} object) a logical indicating whether a column with the strand information have to be added to the result (Defaults to FALSE)
#' 
#' @return
#' A \code{data.frame} with the regions in A. If A was a \code{\link{GRanges}} object, the output will include any metadata present in A.
#' 
#' @seealso \code{\link{toGRanges}}
#' 
#' @examples 
#' A <- data.frame(chr=1, start=c(1, 15, 24), end=c(10, 20, 30), x=c(1,2,3), y=c("a", "b", "c"))
#' 
#' A2 <- toGRanges(A)
#' 
#' toDataframe(A2)
#' 
#' @export toDataframe


#TODO: CHANGE STRANDED TO TRUE OR TO "AUTO", returning a strand column only if there's strand info?

toDataframe <- function(A, stranded=FALSE) {
  
  if(!hasArg(A)) stop("A is missing")
  
  if(is(A, "data.frame")) {
    return(A)
  }
  
  if(is(A, "GRanges")) {
    A <- suppressWarnings(as.data.frame(A)[,-4]) #Return the data in the GRanges removing the strand and width
    names(A)[1] <- "chr"
    if(!stranded) {
      A <- A[,-4]
    }
    return(A)
  }
  
  warning("Unidentified class in toDataFrame. Returning without modifications")
  return(A)
  
}
