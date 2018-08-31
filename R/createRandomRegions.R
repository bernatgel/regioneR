#' Create Random Regions
#' 
#' @description 
#' Creates a set of random regions with a given mean size and standard deviation.
#' 
#' @details
#' A set of nregions will be created and randomly placed over the genome. The lengths of the region set will follow a normal distribution with a mean size \code{length.mean} and a standard deviation \code{length.sd}. The new regions can be made explicitly non overlapping by setting \code{non.overlapping} to TRUE.
#' A mask can be provided so no regions fall in a forbidden part of the genome.
#' 
#' @note
#' If the standard deviation of the length is large with respect to the mean, negative lengths might be created. These region lengths will be
#'  transfromed to into a 1 and so the, for large standard deviations the mean and sd of the lengths are not guaranteed to be the ones in the parameters.
#' 
#' @usage createRandomRegions(nregions=100, length.mean=250, length.sd=20, genome="hg19", mask=NULL, non.overlapping=TRUE)
#' 
#' @param nregions The number of regions to be created.
#' @param length.mean The mean size of the regions created. This is not guaranteed to be the mean of the final region set. See note.
#' @param length.sd The standard deviation of the region size. This is not guaranteed to be the standard deviation of the final region set. See note.
#' @param genome The reference genome to use. A valid genome object. Either a \code{\link{GenomicRanges}} or \code{\link{data.frame}} containing one region per whole chromosome or a character uniquely identifying a genome in \code{\link{BSgenome}} (e.g. "hg19", "mm10" but not "hg"). Internally it uses \code{\link{getGenomeAndMask}}.
#' @param mask The set of regions specifying where a random region can not be (centromeres, repetitive regions, unmappable regions...). A region set in any of the accepted formats (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, ...). \code{\link{NULL}} will try to derive a mask from the genome (currently only works is the genome is a character string) and \code{\link{NA}} explicitly gives an empty mask.
#' @param non.overlapping A boolean stating whether the random regions can overlap (FALSE) or not (TRUE).
#' 
#' @return
#' It returns a \code{\link{GenomicRanges}} object with the regions resulting from the randomization process.
#' 
#' @seealso \code{\link{getGenome}}, \code{\link{getMask}}, \code{\link{getGenomeAndMask}}, \code{\link{characterToBSGenome}}, \code{\link{maskFromBSGenome}}, \code{\link{randomizeRegions}}, \code{\link{resampleRegions}}
#' 
#' @examples
#' genome <- data.frame(c("chr1", "chr2"), c(1, 1), c(180000000, 20000000))
#' mask <- data.frame("chr1", c(20000000, 100000000), c(22000000, 130000000))
#' 
#' createRandomRegions(nregions=10, length.mean=1000, length.sd=500)
#' 
#' createRandomRegions(nregions=10, genome=genome, mask=mask, non.overlapping=TRUE)
#' 
#' @export createRandomRegions


createRandomRegions <- function(nregions=100, length.mean=250, length.sd=20, genome="hg19", mask=NULL, non.overlapping=TRUE) {
  
  if(!is.numeric(nregions)) stop("nregions must be numeric")
  if(!is.numeric(length.mean)) stop("length.mean must be numeric")
  if(!is.numeric(length.sd)) stop("length.sd must be numeric")
  
  gg <- getGenome(genome)
  
  #create a set of regions with the specified length distribution
  lengths<-rnorm(n=nregions, mean=length.mean, sd=length.sd)
  lengths[lengths<1] <- 1
  regs <- data.frame(chr=seqlevels(gg)[1], start=1, end=lengths, stringsAsFactors=FALSE)
  
  return(randomizeRegions(A=regs, genome=genome, mask=mask, allow.overlaps = !non.overlapping))
    
}

