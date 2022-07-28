#' resampleGenome
#'
#' Fast alternative to randomizeRegions. It creates a tiling (binning) of the whole genome 
#' with tiles the mean size of the regions in A and then places the regions by sampling a 
#' length(A) number of tiles and placing the resampled regions there.
#' 
#'
#' @usage resampleGenome(A, simple = FALSE, per.chromosome = FALSE, genome="hg19", min.tile.width=1000, ...)
#'
#' @param A an object of class GenomigRanges
#' @param simple logical, if TRUE the randomization process will not take into account the specific width of each region in A. (defalut = FALSE)
#' @param per.chromosome logical, if TRUE the randomization will be perform by chromosome. (default  = TRUE)
#' @param genome character or GenomicRanges, genome using for the randomization
#' @param min.tile.width integer, the minimum size of the genome tiles. If they are too small, the functions gets very slow and may even fail to work. (default = 1000, 1kb tiles)
#' @param ... further arguments to be passed to other methods.
#'
#'
#' @return  a \code{\link{GenomicRanges}} object. A sample from the \code{universe} with the same length as A.
#'
#' @seealso  \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{randomizeRegions}}, \code{\link{createRandomRegions}}
#'
#'
#' @examples
#'
#' A <- data.frame(chr=1, start=c(2,12,28,35), end=c(5,25,33,43))
#'
#' B <- resampleGenome(A)
#' B
#' width(B)
#' 
#' B2 <- resampleGenome(A, simple=TRUE)
#' B2
#' width(B2)
#'
#' resampleGenome(A, per.chromosome=TRUE)
#'
#'
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges tile
#' @importFrom GenomicRanges resize
#'
#' @export resampleGenome
#'




resampleGenome <- function(A, simple = FALSE, per.chromosome = FALSE, genome = "hg19", min.tile.width=1000, ...) {
  
  if (!methods::hasArg(A)) {
    stop("A is missing")
  }

  if (!is.logical(per.chromosome)) {
    stop("per.chromosome must be logical")
  }
  
  
  A <- toGRanges(A, genome=genome)

  #Build the universe by genome tiling
  mwidth <- round(mean(GenomicRanges::width(A)))
  universe <- unlist(GenomicRanges::tile(getGenome(genome), width = max(mwidth, min.tile.width)))
  
  #Call resample regions
  resampled <- resampleRegions(A=A, universe = universe, per.chromosome = per.chromosome, ...)

  #And resize the selected regions as needed
  if(simple == TRUE) {
    resampled <- GenomicRanges::resize(resampled, width = mwidth, fix = "center", use.names = FALSE)
  } else {
    resampled <- GenomicRanges::resize(resampled, width = GenomicRanges::width(A), fix = "center", use.names = FALSE)
  }

  return(resampled)
}
