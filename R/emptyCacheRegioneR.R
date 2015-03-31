#' Empty Cache regioneR
#'  
#' @description 
#' Empties the caches used by the memoised function in the regioneR package. 
#' 
#' @usage 
#' emptyCacheRegioneR()
#' 
#' @return The cache is emptied
#' 
#' @examples
#' emptyCacheRegioneR()
#' 
#' @export emptyCacheRegioneR


emptyCacheRegioneR <- function() {
  forget(getGenome)
  forget(getMask)
  forget(getGenomeAndMask)
  forget(maskFromBSGenome)
  forget(characterToBSGenome)
  
}
