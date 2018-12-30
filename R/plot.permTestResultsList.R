# Plot Permutation Test Results List
#
# @description
#' Function for plotting the results from a \code{permTestResultsList} object when more than one evaluation function was used.
#' 
#' @method plot permTestResultsList
#' 
#' @param x          an object of class \code{permTestResultsList}.
#' @param ncol       number of plots per row. ncol=NA means ncol=floor(sqrt(length(x)))so the plot is more or less square (default=NA)
#' @param pvalthres  p-value threshold for significance. Default is 0.05.
#' @param plotType   the type of plot to display. This must be one of \code{"Area"} or \code{"Tailed"}. Default is \code{"Area"}.
#' @param main       a character specifying the title of the plot. Defaults to "".
#' @param xlab       a character specifying the label of the x axis. Defaults to NULL, which produces a plot with the evaluation function name as the x axis label.
#' @param ylab       a character specifying the label of the y axis. Defaults to "".
#' @param ...        further arguments to be passed to or from methods.
#' 
#' @return A plot is created on the current graphics device.
#' 
#' @seealso \code{\link{permTest}}
#' 
#' @examples 
#' 
#' genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
#' A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
#' B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
#' 
#' pt <- overlapPermTest(A=A, B=B, ntimes=10, genome=genome, non.overlapping=FALSE)
#' summary(pt)
#' plot(pt)
#' plot(pt, plotType="Tailed")  
#'  
#' pt2 <- permTest(A=A, B=B, ntimes=10, alternative="auto", genome=genome, evaluate.function=list(distance=meanDistance, numberOfOverlaps=numOverlaps), randomize.function=randomizeRegions, non.overlapping=FALSE)
#' summary(pt2)
#' plot(pt2)
#' plot(pt2, plotType="Tailed")
#' 
#' @export 


plot.permTestResultsList<-function(x, ncol=NA, pvalthres=0.05, plotType="Tailed", main="", xlab=NULL, ylab="", ...){
  
  if(!is(x, "permTestResultsList"))  stop("x must be a permTestResultsList object")
  
  if(is.na(ncol)) ncol <- floor(sqrt(length(x)))  
  
  nrow <- ceiling(length(x)/ncol)
  
  old.par <- par(mfrow=c(nrow, ncol))
  
  lapply(x, plot, ...)

  par(mfrow=old.par)
  
}
