#' filterChromosomes
#' listChrTypes
#' 
#' @description 
#'  Prints a list of the available organisms and chromosomes sets in the predefined chromosomes sets information.
#'
#'
#' @usage listChrTypes()
#' 
#' @return the list of available chrs and organisms is printed
#'  
#' @seealso \code{\link{filterChromosomes}}, \code{\link{getChromosomesByOrganism}}
#' 
#' @examples
#' 
#' g <- getGenomeAndMask("hg19")$genome
#' 
#' listChrTypes()
#' 
#' g <- filterChromosomes(g, chr.type="autosomal", organism="hg19") 
#' 
#' @export listChrTypes


listChrTypes <- function() {
  
  chrs <-  getChromosomesByOrganism()
  for(org in names(chrs)) {
    chr.ty <- names(chrs[[org]])
    chr.ty <- chr.ty[chr.ty != "org.name"]
    cat(paste0(chrs[[org]][["org.name"]], " (", org, "): ", paste(chr.ty, collapse=", ")), "\n")
  } 
}

