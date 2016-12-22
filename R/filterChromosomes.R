#' filterChromosomes
#' 
#' @description 
#' Filters the chromosomes in a region set. It can either filter using a predefined chromosome set (e.g. "autosomal
#'  chromosomes in Homo sapiens") or using a custom chromosome set (e.g. only chromosomes "chr22" and "chrX")
#'
#' @usage filterChromosomes(A, organism="hg", chr.type="canonical", keep.chr=NULL)
#' 
#' @param A a region set in any of the formats accepted by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param organism a character indicating the organism from which to get the predefined chromosome sets. It can be the organism code as used in \code{\link{BSgenome}} (e.g. hg for human, mm for mouse...) or the full genome assembly identifier, since any digit will be removed to get the organism code.
#' @param chr.type a character indicating the specific chromosome set to be used. Usually "autosomal" or "canonical", althought other values could be available for certain organisms.
#' @param keep.chr is a character vector stating the names of the chromosomes to keep. Any chromosome not in the vector will be filtered out. If keep.chr is supplied, organism and chr.type are ignored.
#'   
#' @return
#' A \code{\link{GRanges}} object containing only the regions in the original region set belonging to the selected chromosomes. All regions in non selected chromosomes are removed.
#' 
#' @seealso \code{\link{getGenomeAndMask}}, \code{\link{listChrTypes}} \code{\link{getChromosomesByOrganism}}
#' 
#' @examples
#' 
#' g <- getGenomeAndMask("hg19")$genome
#' listChrTypes()
#' g <- filterChromosomes(g, chr.type="autosomal", organism="hg19")
#' g <- filterChromosomes(g, keep.chr=c("chr1", "chr2", "chr3"))
#' 
#' 
#' @export filterChromosomes


filterChromosomes <-function(A,  organism="hg", chr.type="canonical", keep.chr=NULL) { 
  
  A <- toGRanges(A)
  
  if(chr.type == "custom" | !is.null(keep.chr)){
    valid.chr <- keep.chr
  } else {
    org <- getChromosomesByOrganism()
    org.code <- gsub("\\d","", organism) #The name of the organism is assumed to be the assembly identifier minus the digits e.g. hg19 -> hg

    if (org.code %in% names(org)) {
      if (chr.type %in% names(org[[org.code]]) & (chr.type != "org.name")) {
        valid.chr<-as.character(org[[org.code]][chr.type][[1]])
      } else {
        valid.types <- names(org[[org.code]])
        valid.types <- valid.types[valid.types != "org.name"]
        stop(paste("Chromosome type ", chr.type, " for organism", org.code, " (", organism, ") not recognized. Available valid values for ", org.code, " are: ", paste(valid.types, collapse=", ")))
      }
    } else  {
      stop(paste("In filterChromosomes: Organism ", org.code, " (", organism, ") not recognized. Available organisms are: ", paste(names(org), collapse=", ")))
    }
  }
  
  A <- keepSeqlevels(A, valid.chr, pruning.mode="coarse")
  return(A)
  
}