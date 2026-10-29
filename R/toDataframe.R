#' toDataframe
#' 
#' @description
#' Transforms a \code{\link[GenomicRanges]{GRanges}} object or a \code{\link{data.frame}}containing a region set into a \code{\link{data.frame}}.
#' 
#' @details
#' If the oject is of class \code{\link{data.frame}}, it will be returned untouched.
#' 
#' @usage toDataframe(A, stranded=FALSE)
#' 
#' @param A a \code{\link[GenomicRanges]{GRanges}} object.
#' @param stranded    (only used when A is a \code{\link[GenomicRanges]{GRanges}} object) a logical indicating whether a column with the strand information have to be added to the result (Defaults to FALSE)
#' 
#' @return
#' A \code{data.frame} with the regions in A. If A was a \code{\link[GenomicRanges]{GRanges}} object, the output will include any metadata present in A.
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

  if(methods::is(A, "data.frame")) {
    return(A)
  }

  if(methods::is(A, "GRanges")) {
    result <- data.frame(
      chr=as.factor(GenomicRanges::seqnames(A)),
      start=GenomicRanges::start(A),
      end=GenomicRanges::end(A),
      row.names=names(A),
      check.names=FALSE
    )

    if(stranded) {
      result$strand <- as.factor(GenomicRanges::strand(A))
    }

    result <- cbind(
      result,
      as.data.frame(S4Vectors::mcols(A, use.names=FALSE))
    )

    return(result)
  }

  warning("Unidentified class in toDataFrame. Returning without modifications")
  return(A)
  
}
