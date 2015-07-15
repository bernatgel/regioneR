#' Plot a list of localZscore results
#' 
#' @description
#' Function for plotting the a \code{localZScoreResultsList} object.
#' 
#' @method plot localZScoreResultsList
#' 
#' @param lz        an object of class \code{localZScoreResultsList}.
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
#' pt2 <- permTest(A=A, B=B, ntimes=10, randomize.function=randomizeRegions, evaluate.function=list(overlap=numOverlaps, distance=meanDistance), genome=genome, non.overlapping=FALSE)
#' plot(pt2)
#' 
#' lz2 <- localZScore(A=A, B=B, pt2)
#' plot(lz2)
#' 
#' @keywords internal
#' @export


plot.localZScoreResultsList <- function(lz, ncol=NA, main="", num.x.labels=5, ...) {
  
  if(!is(lz, "localZScoreResultsList"))  stop("lz must be a localZScoreResultsList object")
  
  if(is.na(ncol)) ncol <- floor(sqrt(length(lz)))  
  
  nrow <- ceiling(length(lz)/ncol)
  
  old.par <- par(mfrow=c(nrow, ncol))
  
  lapply(lz, plot)
  
  par(mfrow=old.par)
  
}
