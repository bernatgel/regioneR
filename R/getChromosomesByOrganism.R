#' getChromosomesByOrganism
#' 
#' @description 
#' Function to obtain a list of organisms with their canonical and (when applicable) the autosomal chromosome names.
#' This function is not usually used by the end user directly but through the filterChromosomes function.
#' 
#' 
#' @usage getChromosomesByOrganism()
#' 
#' 
#' @return
#' a list with the organism as keys and the list of available chromosome sets as values
#' 
#' @seealso \code{\link{getGenome}}, \code{\link{filterChromosomes}}
#' 
#' @examples
#'
#' chrsByOrg <- getChromosomesByOrganism() 
#' chrsByOrg[["hg"]]
#' chrsByOrg[["hg"]][["autosomal"]]
#' 
#' @export getChromosomesByOrganism
#' @importFrom utils as.roman


getChromosomesByOrganism <-function() {
  
  chromosomesByOrganism <- list(
    hg = list(autosomal=paste0("chr", c(1:22)),
              canonical=paste0("chr",c(1:22,"X","Y")),
              org.name=("Homo sapiens")),
    
    mm = list(autosomal=paste0("chr",c(1:19)),
              canonical=paste0("chr",c(1:19,"X","Y")),
              org.name=("Mus musculus")),
    
    bosTau = list(autosomal=paste0("chr",c(1:29)),
             canonical=paste0("chr",c(1:29,"X","Y")),
             org.name=("Bos taurus")),
    
    ce = list(autosomal=paste0("chr",c(1:5)),
              canonical=paste0("chr",c(1:5,"X","Y")),
              org.name=("Caenorhabditis elegans")),
    
    danRer = list(canonical=paste0("chr",c(1:25)),
                  org.name=("Danio rerio")),
    
    rheMac = list(autosomal=paste0("chr",c(1:20)),
                  canonical=paste0("chr",c(1:20,"X","Y")),
                  org.name=("Macaca mulata")),
    
    rn = list(autosomal=paste0("chr",c(1:20)),
              canonical=paste0("chr",c(1:20,"X","Y")),
              org.name=("Rattus norvegicus")),
    
    sacCer = list(autosomal=paste0("chr",c(as.character(as.roman(1:16)), "M")),
                  canonical=paste0("chr",c(as.character(as.roman(1:16)), "M")),
                  org.name=("Saccharomyces cerevisiae")),
    
    dm = list(autosomal=paste0("chr",c("2L","2R","3L","3R","4")),
              canonical=paste0("chr",c("2L","2R","3L","3R","4","X")),
              org.name=("Drosophila melanogaster")),
    
    panTro = list(autosomal=paste0("chr",c(1,"2A","2B",3:22)),
                  canonical=paste0("chr",c(1,"2A","2B",3:22,"X","Y")),
                  org.name=("Pan troglodytes"))
  )
  
  return(chromosomesByOrganism)
}
