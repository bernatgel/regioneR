#' Plot localZscore results
#' 
#' @description
#' Function for plotting the a \code{localZScoreResults} object.
#' 
#' @method plot localZScoreResults
#' 
#' @param x        an object of class \code{localZScoreResults}.
#' @param main      a character specifying the main title of the plot. Defaults to no title.
#' @param ...       further arguments to be passed to or from methods.
#' 
#' @return A plot is created on the current graphics device.
#' 
#' @seealso \code{\link{localZScore}}
#' 
#' @examples 
#'
#' genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
#' A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
#' B <- c(A, createRandomRegions(nregions=10, length.mean=100000, length.sd=20000, genome=genome, non.overlapping=FALSE))
#' 
#' pt <- overlapPermTest(A=A, B=B, ntimes=10, genome=genome, non.overlapping=FALSE)
#'  
#' lz <- localZScore(A=A, B=B, pt=pt)
#' plot(lz)
#' 
#' @keywords internal
#' @export


plot.localZScoreResults <- function(x, main="", ...) {
  toLabel <- function(n) {
    if(abs(n) < 1000) return(as.character(n))
    if(abs(n) < 1000000) return(paste0(as.character(round(n/10)/100), "Kb"))
    return(paste0(as.character(round(n/10000)/100), "Mb"))
  }

  if(nchar(main)==0) main <- "Local z-score"


  old.scipen <- options("scipen")
  options(scipen=999)
  x.labs <- sapply(x$shifts, toLabel)
  y.max <- max(x$shifted.z.scores, 2)
  y.min <- min(x$shifted.z.scores, -2)
  plot(x=x$shifts, y=x$shifted.z.scores, type="l", ylim=c(y.min, y.max), ylab="Shifted z-scores", xlab="Shifts", main=main, xaxt="n", las=1, ...)
  axis(1, at=x$shifts, labels=x.labs, las=2, cex.axis=0.7, tck=-.01, ...)
  box(lwd=1.2)
  options(scipen=old.scipen)
  
}
