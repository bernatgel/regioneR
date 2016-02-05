#' characterToBSGenome
#' 
#' @description 
#' Given a character string with the "name" of a genome, it returns a \code{\link{BSgenome}} object if available.
#' 
#' @note
#' This function is memoised (cached) using the \code{memoise} package. To empty the cache, use \code{forget(charecterToBSGenome)}
#' 
#'  @usage characterToBSGenome(genome.name) 
# @usage characterToBSGenome(...)
#' 
#' @param ... a genome.name parameter is needed, a character string uniquely identifying a \code{\link{BSgenome}} (e.g. "hg19", "mm10" are ok, but "hg" is not)
#' 
#' @return
#' A \code{\link{BSgenome}} object
#' 
#' @examples
#' g <- characterToBSGenome("hg19")

#' @seealso \code{\link{getGenomeAndMask}}, \code{\link{maskFromBSGenome}}
#' 
#' @export characterToBSGenome


characterToBSGenome <- memoise(function(genome.name) {
  
  if(!hasArg(genome.name)) stop("parameter genome.name is required")
  if(!is.character(genome.name)) stop("genome.name must be a character")

  bsg <- NULL
  #Try to get the masked BSgenome with the getBSgenome
  tryCatch(
    expr={
      bsg <- getBSgenome(genome.name, masked=TRUE)
    },
    error = function(err) {
     #do nothing 
    })
  
  
  if(is.null(bsg)) { #Try to get the unmasked BSgenome with the getBSgenome if the masked was not available
    bsg <- getBSgenome(genome.name, masked=FALSE)
    warning(paste0("The masked version of '", genome.name, "' is not installed. Using the unmasked version. This means that no automatic masking will be available."))
  }
  
  return(bsg)
})


