#' Empty Cache regioneR
#'  
#' @description 
#' Empties the caches used by the memoised functions in the regioneR package. 
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
#' 
#' @importFrom memoise forget


emptyCacheRegioneR <- function() {
  memoise::forget(getGenome)
  memoise::forget(getMask)
  memoise::forget(getGenomeAndMask)
  memoise::forget(maskFromBSGenome)
  memoise::forget(characterToBSGenome)
  
}
