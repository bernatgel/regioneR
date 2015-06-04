#' Plot localZscore results
#' 
#' @description
#' Function for plotting the a \code{localZScoreResults} object.
#' 
#' @method plot localZScoreResults
#' 
#' @param lz        an object of class \code{localZScoreResults}.
#' @param main      a character specifying the main title of the plot. Defaults to no title.
#' @param num.x.labels  a numeric specifying the number of ticks to label the x axis. The total number will be 2*num.x.labels + 1. Defaults to 5.
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


plot.localZScoreResults <- function(lz, main="", num.x.labels=5, ...) {
  #Convert a number to a "human readable" label
  toLabel <- function(n) {
    if(abs(n) < 1000) return(as.character(n))
    if(abs(n) < 1000000) return(paste0(as.character(round(n/10)/100), "Kb"))
    return(paste0(as.character(round(n/10000)/100), "Mb"))
  }

  if(nchar(main)==0) main <- "Local z-score"

  old.scipen <- options("scipen")
  options(scipen=999)
  
  #Set the positions for the x labels
  if(num.x.labels < 1) {
    x.lab.pos <- 0
  } else {
    x.lab.dist <- floor(lz$window/num.x.labels)
    x.lab.pos <- (1:num.x.labels)*x.lab.dist
    x.lab.pos <- c(rev(-1*x.lab.pos), 0, x.lab.pos)
  }
  x.labs <- sapply(x.lab.pos, toLabel)
  
  y.max <- max(lz$shifted.z.scores, 2)
  y.min <- min(lz$shifted.z.scores, -2)
  plot(x=lz$shifts, y=lz$shifted.z.scores, type="l", ylim=c(y.min, y.max), ylab="Shifted z-scores", xlab="Shifts", main=main, xaxt="n", las=1, ...)
  if(num.x.labels != 0) {
    axis(1, at=x.lab.pos, labels=x.labs, las=2, cex.axis=0.7, tck=-.01, ...)
  }
  box(lwd=1.2)
  
  options(scipen=old.scipen)
  
}
