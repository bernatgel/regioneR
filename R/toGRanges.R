#' toGRanges
#' 
#' @description
#' Transforms a file or an object containing a region set into a \code{\link{GRanges}} object. 
#' 
#' @details
#' If A is already a \code{\link{GRanges}} object, it will be returned untouched. 
#' 
#' If A is a file name or connection to a file in any of the formats supported by \code{rtracklayer}'s import function (BED, GFF...) 
#' it will be imported using \code{rtracklayer}.
#'  
#' If A is a data frame, the function will assume the first three columns are chromosome, start and end and create a \code{\link{GRanges}} object. Any additional 
#' column will be considered metadata and stored as such in the \code{\link{GRanges}} object. 
#' 
#' If A is not a data.frame and there are more parameters, it will try to build a data.frame with all 
#' parameters and use that data.frame to build the GRanges. This allows the user to call it like
#'  \code{toGRanges("chr1", 10, 20)}.
#'  
#' If A is a character or a character vector and it's not a file or a URL, it 
#' assumes it's a genomic position description in the form used by UCSC or 
#' IGV, "chr2:1000-2000". It will try to parse the character strings into
#' chromosome, start and end and create a GRanges. The parser can deal with 
#' commas separating thousands (e.g. "chr2:1,000-2,000") and with the comma 
#' used as a start/end separator (e.g. "chr2:1000,2000"). These different 
#' variants can be mixed in the same character vector.
#' 
#' The \code{genome} parameter can be used to set the genome information of
#' the created GRanges. It can be either a \code{\link{BSgenome}} object or a 
#' character string defining a genome (e.g. "hg19", "mm10"...) as accepted 
#' by the \code{BSgenome::getBSgenome} function. If a valid genome is
#' given and the corresponding BSgenome package is installed, the genome
#' information will be attached to the GRanges. If the chromosome naming style
#' from the GRanges and the genome object are different, it will try to change 
#' the GRanges styles to match those of the genome using 
#' \code{GenomeInfoDb::seqlevelsStyle}.
#' 
#' 
#' 
#' @usage toGRanges(A, ..., genome=NULL)
#' 
#' @param A  a \code{\link{data.frame}} containing a region set, a \code{\link{GRanges}} object, a BED file or any type of file supported by \code{rtracklayer}. If there are more than 1 argument, it will build a dataframe out ouf them and process it as usual. If there's only a single argument and it's a character, if it's not an existing file name it will be treated as the definition of a genomic region in the UCSC/IGV format (i.e. "chr9:34229289-34982376") and parsed. 
#' @param ... further arguments to be passed to other methods. 
#' @param genome (character or BSgenome) The genome info to be attached to the created GRanges. If NULL no genome info will be attached. (defaults to NULL)
#' 
#' @return
#' A \code{\link{GRanges}} object with the regions in A
#' 
#' @seealso \code{\link{toDataframe}}
#' 
#' @examples
#' A <- data.frame(chr=1, start=c(1, 15, 24), end=c(10, 20, 30),  x=c(1,2,3), y=c("a", "b", "c"))
#' gr1 <- toGRanges(A)
#' 
#' #No need to give the data.frame columns any specific name
#' A <- data.frame(1, c(1, 15, 24), c(10, 20, 30),  x=c(1,2,3), y=c("a", "b", "c"))
#' gr2 <- toGRanges(A)
#' 
#' #We can pass the data without building the data.frame
#' gr3 <- toGRanges("chr9", 34229289, 34982376, x="X")
#' 
#' #And each argument can be a vector (they will be recycled as needed)
#' gr4 <- toGRanges("chr9", c(34229289, 40000000), c(34982376, 50000000), x="X", y=c("a", "b"))
#' 
#' #toGRanges will automatically convert the second and third argument into numerics
#' gr5 <- toGRanges("chr9", "34229289", "34982376") 
#' 
#' #It can be a file from disk
#' bed.file <- system.file("extdata", "my.special.genes.bed", package="regioneR")
#' gr6 <- toGRanges(bed.file)
#' 
#' #Or a URL to a valid file
#' gr7 <- toGRanges("http://molb7621.github.io/workshop/_downloads/lamina.bed")
#' 
#' #It can also parse genomic location strings
#' gr8 <- toGRanges("chr9:34229289-34982376")
#' 
#' #more than one
#' gr9 <- toGRanges(c("chr9:34229289-34982376", "chr10:1000-2000"))
#' 
#' #even with mixed strange and mixed syntaxes
#' gr10 <- toGRanges(c("chr4:3873-92928", "chr4:3873,92928", "chr5:33,444-45,555"))
#' 
#' #if the genome is given it is used to annotate the resulting GRanges
#' gr11 <- toGRanges(c("chr9:34229289-34982376", "chr10:1000-2000"), genome="hg19")
#' 
#' 
#' #and the genome is added to the GRanges even if A is a GRanges
#' gr12 <- toGRanges(gr6, genome="hg19")
#' 
#' #And it will change the chromosome naming of the GRanges to match that of the
#' #genome if it is possible (using GenomeInfoDb::seqlevelsStyle)
#' gr2
#' gr13 <- toGRanges(gr2, genome="hg19")
#' 
#' 
#' 
#' @export toGRanges
#' 
#' @importFrom GenomicRanges GRanges elementMetadata
#' @import BSgenome
#' @importFrom rtracklayer import
#' @importFrom utils read.delim read.csv head
#' @import parallel
#' @import GenomeInfoDb
#' 
#' 
#' 


toGRanges <- function(A, ..., genome=NULL) {
    
  if(!hasArg(A)) stop("A is missing")
  
  if(methods::is(A, "GRanges")) {
    return(setGenomeToGRanges(A, genome))
  }
  
  
  #if there are more than one parameters, try to build a data.frame from them
  if(length(list(...))>0) {
    tryCatch(A <- do.call(data.frame, c(list(A), list(...), list(stringsAsFactors=FALSE))),
             error=function(e){}, warning=function(e){})
  } else { #If it's single parameter
    #If a character, it's not a file on disk and does not contain "://" (so it's not a URL), try to parse it as a genome coordinates string
    if(all(is.character(A)) && !any(file.exists(A)) && !any(grepl(pattern = "://", A))) {
      #If there are no "-" and only one "," may be in this format "chr:start,end". Replace the "," with a "-" and proceed.
      A2 <- A
      comma.sep <- !grepl("-", A2) & (nchar(A2) - nchar(gsub(",", "", A2))==1)
      A2[comma.sep] <- gsub(",", "-", A2[comma.sep])
      #Remove the "," used as thousand separators
      A2 <- gsub(pattern = ",", replacement = "", x = A2)
      
      #Split on : and -
      A2 <- strsplit(A2, "[:-]+")
      #And count them
      if(!all(lapply(A2, length)==3)) {
        failing <- A[which(lapply(A2, length)!=3)]
        stop("Error when parsing the genomic region definition strings. There are ", length(failing), " malformed strings: ", paste0(paste0(head(failing, n=5), collapse = ","), ifelse(length(failing)>5, paste0(" and ", length(failing)-5, " more."), ".")))
      } else {
        A <- data.frame(do.call(rbind, A2))
        A[,2] <- as.numeric(as.character(A[,2])) #we need the as.character because they are factors
        A[,3] <- as.numeric(as.character(A[,3]))
      }
    }
  }
  
  #if it's a dataframe, assume the three first columns are: chr, start and end
  if(methods::is(A, "data.frame")) {
    if(length(A)==0) stop("In GRanges: A cannot be a data.frame with 0 columns.")
      
    if(length(A)==1) { #If the data.frame has inly one column, treat it as a character vector
      return(toGRanges(A[,1]), genome=genome)
    }
    
    #If the data.frame has only two columns, repeat the second column to create "one-base-wide" regions
    if(length(A)==2) {
      A <- cbind(A, A[,2], stringsAsFactors=FALSE)
    }
    
    chrs <- as.character(A[,1]) #Transform the first column into a character. It does not work with factors.
    #and transform the second and third into numerics
    if(is.factor(A[,2])) A[,2] <- as.character(A[,2])
    if(is.factor(A[,3])) A[,3] <- as.character(A[,3])
    start <- as.numeric(A[,2])
    end <- as.numeric(A[,3])
    
    gr <- GenomicRanges::GRanges(seqnames=chrs, ranges=IRanges::IRanges(start=start, end=end))
    #We cannot assign the metadata in a single line because when only one metadata column was requested, 
    #it was automatically transformed into a vector and it lost its name. 
    if(ncol(A)>3) { #if there's metadata or strand information
      #Store the metadata
      metadata.columns <- c(4:ncol(A))
      if(length(metadata.columns>0)) {
        original.metadata <- as.data.frame(A[,metadata.columns])   
        names(original.metadata) <- names(A)[metadata.columns]
        GenomicRanges::elementMetadata(gr) <- original.metadata
      }
    }
  } else { #It's not a data.frame

    #If its something else, try with rtracklayer and see if it can read and import from it.
    #If it fails and raises an error (usually because unknown format (txt)), try to open it with read.delim and if it fails, try with read.csv
    gr <- tryCatch(
      expr = {rtracklayer::import(con=A, ...)},
      error = function(err) {
        return(tryCatch(
          expr = {toGRanges(utils::read.delim(A, ...))},
          error = function(err) {
            return(tryCatch(
              expr = {toGRanges(utils::read.csv(A, ...))},
              error = function(e) {stop("Error in toGRanges when trying to read the file \"", A, "\": ", e)}
            ))
          }
        ))
      }
    )
  }
  
  return(setGenomeToGRanges(gr, genome))
}  
  
  

setGenomeToGRanges <- function(gr, genome) {
  if(!is.null(genome)) {
    if(is.character(genome)) genome <- characterToBSGenome(genome)
    if(!methods::is(genome, "BSgenome")) {
      warning("Invalid 'genome' argument. Ignoring.")
      genome <- NULL
    }
  }
  if(!is.null(genome)) {
    GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(genome)
    GenomeInfoDb::seqlevels(gr, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(genome)
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(genome)
  }
  return(gr)
}

