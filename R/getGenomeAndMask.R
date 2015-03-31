#' getGenomeAndMask
#' 
#' @description 
#' Function to obtain a valid genome and mask pair given a valid genome identifier and optionally a mask. 
#' 
#' If the genome is not a \code{\link{BSgenome}} object or a character string uniquely identifying a \code{\link{BSgenome}} package installed, it will return the genome "as is". If a mask is provided, it will simply return it. Otherwise it will return the mask returned by \code{\link{getMask}(genome)} or an empty mask if genome is not a valid \code{\link{BSgenome}} or \code{\link{BSgenome}} identifier.
#' 
#' @note
#' This function is memoised (cached) using the \code{\link{memoise}} package. To empty the cache, use \code{\link{forget}(getGenomeAndMask)}
#' 
#' 
# @usage getGenomeAndMask(genome, mask=NULL) <- Real Documentation. Problems with codoc
#' @usage getGenomeAndMask(...)
#' 
#' @param ...  a genome parameter is required, the genome object or genome identifier. A mask parameter is optional: the mask of the genome. If mask is \code{\link{NULL}}, it will try to get a mask from the genome. If mask is \code{\link{NA}} it will return an empty mask.
#' 
#' @return
#' A list with two elements: genome and mask. Genome and mask are GRanges objects. 
#' 
#' @seealso \code{\link{getMask}}, \code{\link{getGenome}}, \code{\link{characterToBSGenome}}, \code{\link{maskFromBSGenome}}, \code{\link{emptyCacheRegioneR}}
#' 
#' @examples
#' getGenomeAndMask("hg19", mask=NA)
#' 
#' getGenomeAndMask(genome=data.frame(c("chrA", "chrB"), c(15000000, 10000000)), mask=NA)
#' 
#' @export getGenomeAndMask

getGenomeAndMask <- memoise(function(genome, mask=NULL) {

  
  #if genome is a character, get it from the BS packages
  if(is.character(genome)) { 
    genome <- characterToBSGenome(genome)
  }
  
  
  if(is(genome, "MaskedBSgenome") && is.null(mask)) {
    mask <- getMask(genome)
  } else {
    #check if it seems to be a valid mask
      # try to create a GRanges object with it, if it works, assume its valid. If not, return an empty mask
      mask <- try(exp=toGRanges(mask), silent=TRUE)
      if(!(is(mask, "GenomicRanges"))) {mask <- toGRanges(data.frame(chr=character(), start=numeric(), end=numeric())) }
      
  }
  
  genome <- getGenome(genome)
  
  return(list(mask=mask, genome=genome))
})
