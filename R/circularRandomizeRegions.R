#' Circular Randomize Regions
#' 
#' @description 
#' Given a set of regions A and a genome, this function returns a new set of regions created by applying a random
#' spin to each chromosome. 
#' 
#' @details
#' This randomization strategy is useful when the spatial relation between the regions in the RS is important and has to be conserved.
#' 
#' @usage 
#' circularRandomizeRegions(A, genome="hg19", mask=NULL, max.mask.overlap=NULL, max.retries=10, verbose=TRUE, ...)
#' 
#' @param A The set of regions to randomize. A region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param genome The reference genome to use. A valid genome object. Either a \code{\link{GenomicRanges}} or \code{\link{data.frame}} containing one region per whole chromosome or a character uniquely identifying a genome in \code{\link{BSgenome}} (e.g. "hg19", "mm10" but not "hg"). Internally it uses \code{\link{getGenomeAndMask}}.
#' @param mask The set of regions specifying where a random region can not be (centromeres, repetitive regions, unmappable regions...). A region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}},\code{\link{data.frame}}, ...). If \code{\link{NULL}} it will try to derive a mask from the genome (currently only works is the genome is a character string) and if \code{\link{NA}} it will explicitly give an empty mask.
#' @param max.mask.overlap numeric value
#' @param max.retries numeric value
#' @param verbose a boolean.
#' @param ... further arguments to be passed to or from methods.
#' 
#' 
#' @return
#' It returns a \code{\link{GenomicRanges}} object with the regions resulting from the randomization process.
#' 
#' @seealso  \code{\link{randomizeRegions}}, \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{getGenome}}, \code{\link{getMask}}, \code{\link{getGenomeAndMask}}, \code{\link{characterToBSGenome}}, \code{\link{maskFromBSGenome}}, \code{\link{resampleRegions}}, \code{\link{createRandomRegions}}
#' 
#' @examples
#' A <- data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40))
#' 
#' mask <- data.frame("chr1", c(20000000, 100000000), c(22000000, 130000000))
#' 
#' genome <- data.frame(c("chr1", "chr2"), c(1, 1), c(180000000, 20000000))
#' 
#' circularRandomizeRegions(A)
#' 
#' circularRandomizeRegions(A, genome=genome, mask=mask, per.chromosome=TRUE, non.overlapping=TRUE)
#' 
#' @export circularRandomizeRegions
#' 




circularRandomizeRegions <- function(A, genome="hg19", mask=NULL, max.mask.overlap=NULL, max.retries=10, verbose=TRUE, ...) {
  A <- toDataframe(toGRanges(A))
  gam <- getGenomeAndMask(genome, mask)
  
  
  
  getRandomChr <- function(chr, gam) {
	
    chr.A <- A[A$chr==chr,]
    chr.len <- end(gam$genome[seqnames(gam$genome)==as.character(chr)])
	
    spin <- floor(runif(1)*chr.len)
    spin.A <- spinChromosome(chr.A, spin, chr.len)
    if(!is.null(max.mask.overlap)) {
      #check mask
      num.ov <- numOverlaps(A=spin.A, B=gam$mask, count.once=TRUE)

      num.retries <- 0
      while(num.ov/nrow(spin.A) > max.mask.overlap & num.retries < max.retries) {
        if(verbose==TRUE) {
          message(paste0("Chromosome ", chr, ". Too much overlap with the mask: Num. Regions: ", nrow(chr.A), "  Num. Overlaps: ", num.ov, "  Pct. Overlap: ", (num.ov/nrow(spin.A))*100, ". Retrying..."))
        }
        spin <- floor(runif(1)*chr.len)
        spin.A <- spinChromosome(chr.A, spin, chr.len)
        num.ov <- numOverlaps(A=spin.A, B=gam$mask, count.once=TRUE)
        num.retries <- num.retries + 1
      }
      if(num.retries >= max.retries) {
        stop(paste0("ERROR: After ", max.retries, " retries, it was not possible to find a spin for chromosome ", chr," with an overlap with the mask lower than ", max.mask.overlap, ". The mask is too dense or there are too many features."))
      }
    }
    return(spin.A) 
  }
  
  rand.A <- lapply(unique(A$chr), getRandomChr, gam=gam)
  rand.A <- do.call(rbind, rand.A)
  return(toGRanges(rand.A))
  
}




spinChromosome <- function(A, spin, chr.len) {

  A$start <- (A$start + spin) %% chr.len
  A$end <- (A$end + spin) %% chr.len
  
  part.out <- A$end < A$start
  if(any(part.out)) {
    C <- A[part.out,]
    A[part.out, "end"] <- chr.len
    C$start <- rep(1, nrow(C)) 
  
    A <- rbind(A, C)
  }

  return(A)
}



