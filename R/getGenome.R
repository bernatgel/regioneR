#' getGenome
#' 
#' @description 
#' 
#' Function to obtain a genome
#' 
#' @return
#' A GRanges object with the "genome" data c(Chromosome, Start (by default, 1), Chromosome Length) given a \code{\link{BSgenome}}, a genome name, a \code{\link{data.frame}} or a GRanges.
#' 
#' @details 
#'  
#' If genome is a \code{\link{BSgenome}} (from the package \code{BioStrings}), it will transform it into a \code{\link{GRanges}} with chromosomes and chromosome lengths.
#' 
#' If genome is a \code{\link{data.frame}} with 3 columns, it will transform it into a GRanges.
#' 
#' If genome is a \code{\link{data.frame}} with 2 columns, it will assume the first is the chromosome, the second is the length of the chromosomes and will add 1 as start.
#' 
#' If genome is a \code{character} string uniquely identifying a \code{\link{BSgenome}} installed in the system (e.g. "hg19", "mm10",... but not "hg"), it will create a genome based on the \code{\link{BSgenome}} object identified by the character string.
#' 
#' If genome is a \code{ \link{GRanges}} object, it will return it as is.
#' 
#' If genome is non of the above, it will give a warning and try to transform it into a GRanges using \link{toGRanges}. This can be helpful if \code{genome} is a connection to a file.
#' 
# @usage getGenome(genome) <- Real Documentation. Problems with codoc
#' 
#' @usage getGenome(genome)
#' 
#' @param genome The genome object or genome identifier.
#' 
#' @return
#' A \code{\link{GRanges}} representing the genome with one region per chromosome.
#' 
#' @note 
#' 
#' This function is memoised (cached) using the \code{\link{memoise}} package. To empty the cache, use \code{\link{forget}(getGenome)}
#' 
#' Please note that passing this function the path to a file will not work, since it will assume the character is the identifier of a genome. To read the genome
#' from a file, please use \code{getGenome(toGRanges("path/to/file"))}
#' 
#' @seealso \code{\link{getMask}}, \code{\link{getGenomeAndMask}}, \code{\link{characterToBSGenome}}, \code{\link{maskFromBSGenome}}, \code{\link{emptyCacheRegioneR}}
#' 
#' @examples
#' getGenome("hg19")
#' 
#' getGenome(data.frame(c("chrA", "chrB"), c(15000000, 10000000)))
#'  
#' @export getGenome
#' 


getGenome <- memoise(function(genome) {

    if(!hasArg(genome)) {stop("No genome was specified. genome is a required parameter")}
  
    if(is(genome, "GRanges")) {
      return(genome)
    }
    
    #if specified as a character, get it from the BS packages
    if(is.character(genome)) { 
      genome <- characterToBSGenome(genome)
    }
    
    if(is(genome, "BSgenome")) { #it may be a BS genome because it was originally or because it has been transformed from a chracter
      ss <- seqinfo(genome)
      genome <- data.frame(chr=as.character(ss@seqnames), start=1, end=as.numeric(ss@seqlengths), stringsAsFactors=FALSE)
      
      return(toGRanges(genome))
    }
     
      
    #if the genome is a data frame (not GRanges) and has no starts but only lengths, add them
    if(is(genome, "data.frame") && dim(genome)[2]==2) { 
      genome <- data.frame(chr=genome[,1], chr.start=rep(1,length(genome[,1])), chr.end=as.numeric(genome[,2]), stringsAsFactors=FALSE)
     
      return(toGRanges(genome))
    }
  
    if(is(genome, "data.frame") && dim(genome)[2]==3) {
      return(toGRanges(genome))
    }
    
    warning("Genome format not identified. Trying to to transform with toGRanges.")
    return(toGRanges(genome))
})
