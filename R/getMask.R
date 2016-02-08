#' getMask
#' 
#' @description 
#' Function to obtain a mask given a genome available as a \code{\link{BSgenome}}. The mask returned is the merge of all the active masks in the \code{\link{BSgenome}}.
#' 
#' Since it uses \code{\link{characterToBSGenome}}, the genome can be either a \code{\link{BSgenome}} object or a character string uniquely identifying the a \code{\link{BSgenome}} object installed.
#' 
#' @note
#' This function is memoised (cached) using the \code{\link{memoise}} package. To empty the cache, use \code{\link{forget}(getMask)}
#' 
#' @usage getMask(genome)
# @usage getMask(...)
#' 
#' @param genome The genome from where the mask will be extracted. It can be either a \code{\link{BSgenome}} object or a character string uniquely identifying a \code{\link{BSgenome}} object installed (e.g. "hg19", "mm10", ...)
#' 
#' @return
#' A \code{\link{GRanges}} object with the genomic regions to be masked out
#' 
#' @seealso \code{\link{getGenome}}, \code{\link{getGenomeAndMask}}, \code{\link{characterToBSGenome}}, \code{\link{maskFromBSGenome}}, \code{\link{emptyCacheRegioneR}}
#' 
#' @examples
#' hg19.mask <- getMask("hg19")
#' 
#' hg19.mask
#' 
#' @export getMask


getMask <- memoise(function(genome) {
 
  mask <- NULL
  #if specified as a character, get it from the BS packages
  if(is.character(genome)) { 
    genome <- characterToBSGenome(genome)
  }
  
  if(is(genome, "BSgenome")) { #it may be a BS genome because it was originally or because it has been transformed from a chracter
    mask <- maskFromBSGenome(genome)
  }
  
  return(mask)
})
