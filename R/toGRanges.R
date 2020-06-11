#' toGRanges
#' 
#' @description
#' Transforms a file or an object containing a region set into a 
#' \code{\link{GRanges}} object. 
#' 
#' @details
#' If A is already a \code{\link{GRanges}} object, it will be returned untouched. 
#' 
#' If A is a data frame, the function will assume the first three columns are 
#' chromosome, start and end and create a \code{\link{GRanges}} object. Any 
#' additional column will be considered metadata and stored as such in the 
#' \code{\link{GRanges}} object. There are 2 special cases: 1) if A is a 
#' data.frame with only 2 columns, it will assume the first one is the
#' chromosome and the second one the position and it will create a GRanges with 
#' single base regions and 2) if the data.frame has the first 3 columns named
#' "SNP", "CHR" and "BP" it will shuffle the columns and repeat "BP" to build
#' a GRanges of single base regions (this is the standard ouput format of plink).
#' 
#' If A is not a data.frame and there are more parameters, it will try to build 
#' a data.frame with all parameters and use that data.frame to build the 
#' GRanges. This allows the user to call it like 
#' \code{toGRanges("chr1", 10, 20)}. 
#'  
#' If A is a character or a character vector and it's not a file or a URL, it 
#' assumes it's a genomic position description in the form used by UCSC or 
#' IGV, "chr2:1000-2000". It will try to parse the character strings into
#' chromosome, start and end and create a GRanges. The parser can deal with 
#' commas separating thousands (e.g. "chr2:1,000-2,000") and with the comma 
#' used as a start/end separator (e.g. "chr2:1000,2000"). These different 
#' variants can be mixed in the same character vector.
#' 
#' If A is a "SimpleRleList" it will be interpreted as the result from 
#' GenomicRanges::coverage and the function will return a GRanges with a
#' single metadata column named "coverage".
#' 
#' If A is a file name (local or remote) or a connection to a file, it will try
#' to load it in different ways:
#'   * BED files (identified by a "bed" extension): will be loaded using 
#'   \code{rtracklayer::import} function. Coordinates are 0 based as
#'   described in the BED specification (https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
#'   * PLINK assoc files (identified by ".assoc", ".assoc.fisher", 
#'   ".assoc.dosage", ".assoc.linear", ".assoc.logistic"): will be loaded 
#'   as single-base ranges with all original columns present and the SNPs ids
#'   as the ranges names
#'   * Any other file: It assumes the file is a "generic" tabular file. To load
#'    it it will ignore any header line starting with \code{comment.char}, 
#'    autodetect the field separator (if not provided by the user), 
#'    autodetect if it has a header and read it accordingly.
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
#' @note **IMPORTANT:** Regarding the coordinates, BED files are 0 based 
#' while \code{data.frames} and generic files are treated as 1 based. Therefore
#' reading a line "chr9  100   200" from a BED file will create a 99 bases wide
#' interval starting at base 101 and ending at 200 but reading it from a txt 
#' file or from a \code{data.frame} will create a 100 bases wide interval 
#' starting at 100 and ending at 200. This is specially relevant in 1bp
#' intervals. For example, the 10th base of chromosome 1 would be 
#' "chr1  9  10" in a BED file and "chr1  10  10" in a txt file.
#' 
#' @usage toGRanges(A, ..., genome=NULL, sep=NULL, comment.char="#")
#' 
#' @param A  a \code{\link{data.frame}} containing a region set, a \code{\link{GRanges}} object, a BED file, any type of file supported by \code{rtracklayer::import} or a \code{"SimpleRleList"} returned by \code{GenomicRanges::coverage}. If there are more than 1 argument, it will build a dataframe out ouf them and process it as usual. If there's only a single argument and it's a character, if it's not an existing file name it will be treated as the definition of a genomic region in the UCSC/IGV format (i.e. "chr9:34229289-34982376") and parsed. 
#' @param ... further arguments to be passed to other methods. 
#' @param genome (character or BSgenome) The genome info to be attached to the created GRanges. If NULL no genome info will be attached. (defaults to NULL)
#' @param sep (character) The field separator in the text file. If NULL it will be automatically guessed. Only used when reading some file formats. (Defaults to NULL)
#' @param comment.char (character) The character marking comment lines. Only used when reading some file formats. (Defaults to "#")
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
#' bed.file <- system.file("extdata", "my.special.genes.txt", package="regioneR")
#' gr6 <- toGRanges(bed.file)
#' 
#' #Or a URL to a valid file
#' #gr7 <- toGRanges("http://path.to/my.bed")
#' 
#' #It can also parse genomic location strings
#' gr8 <- toGRanges("chr9:34229289-34982376")
#' 
#' #more than one
#' gr9 <- toGRanges(c("chr9:34229289-34982376", "chr10:1000-2000"))
#' 
#' #even with strange and mixed syntaxes
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
#' #in addition, it can convert other objects into GRanges such as the 
#' #result of GenomicRanges::coverage
#' 
#' gr14 <- toGRanges(c("1:1-20", "1:5-25", "1:18-40"))
#' cover <- GenomicRanges::coverage(gr14)
#' gr15 <- toGRanges(cover)
#' 
#' 
#' 
#' @export toGRanges
#' 
#' @importFrom GenomicRanges GRanges elementMetadata GRangesList
#' @import BSgenome
#' @importFrom rtracklayer import
#' @importFrom tools file_ext
#' @importFrom IRanges tolower
#' @importFrom utils head
#' @import parallel
#' @import GenomeInfoDb
#' 
#' 
#' 


toGRanges <- function(A, ..., genome=NULL, sep=NULL, comment.char="#") {
    
  if(!hasArg(A)) stop("A is missing")
  
  #If it's a GRanges, set the correct genome if needed and return
  if(methods::is(A, "GRanges")) {
    return(setGenomeToGRanges(A, genome))
  }
  
  #If it's a SimpleRleList, assume it's a coverage object as the one generated
  #by GenomicRanges::coverage
  if(methods::is(A, "SimpleRleList")) {
    return(setGenomeToGRanges(coverageToGRanges(A), genome))
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
      
    if(length(A)==1) { #If the data.frame has only one column, treat it as a character vector
      return(toGRanges(A[,1]), genome=genome)
    }
    
    #If the data.frame has only two columns, repeat the second column to create "one-base-wide" regions
    if(length(A)==2) {
      A <- cbind(A, A[,2], stringsAsFactors=FALSE)
    }
    
    #Assoc files fom plink and the like: they have specific column names
    if(ncol(A)>=3 && all(names(A)[1:3] %in% c("CHR", "BP", "SNP"))) {
      if(ncol(A)>3) {
        A <- A[, c("CHR", "BP", "BP", "SNP", names(A)[4:ncol(A)])]
      } else {
        A <- A[, c("CHR", "BP", "BP", "SNP")]        
      }
      if(!any(duplicated(A[,4]))) {
        row.names(A) <- A[,4]
      }
    }
    
    chrs <- as.character(A[,1]) #Transform the first column into a character. It does not work with factors.
    #and transform the second and third into numerics
    if(is.factor(A[,2])) A[,2] <- as.character(A[,2])
    if(is.factor(A[,3])) A[,3] <- as.character(A[,3])
    start <- as.numeric(A[,2])
    end <- as.numeric(A[,3])
    
    gr <- GenomicRanges::GRanges(seqnames=chrs, ranges=IRanges::IRanges(start=start, end=end))
    
    
    #assign the rownames of A to the GRanges too (if not duplicated!)
    if(length(gr)>0 && !any(duplicated(row.names(A)))) {
      tryCatch(names(gr) <- rownames(A), error=function(e){})
    }
    
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
    #If we are here, assume it's a file
    gr <- fileToGRanges(A, ..., sep=sep, genome=genome)
  }
  
  return(setGenomeToGRanges(gr, genome))
}  
  
  

setGenomeToGRanges <- function(gr, genome) {
  if(!is.null(genome) && is.character(genome)) genome <- characterToBSGenome(genome[1])
  
  if(!is.null(genome)) {
    if(methods::is(genome, "BSgenome")) {
      #Seems weird, but at least the  "BSgenome.Hsapiens.1000genomes.hs37d5" returns two styles. Use only the first one
      GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(genome)[1] 
      GenomeInfoDb::seqlevels(gr, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(genome)
      GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(genome)  
    }
    if(methods::is(genome, "SeqInfo")) {
      GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(genome)[1] 
      GenomeInfoDb::seqlevels(gr, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(genome)
      GenomeInfoDb::seqinfo(gr) <- genome
    }
  }
  return(gr)
}

#Transform the SimpleRleList generated by GenomicRanges::coverage into a 
#Granges with the coverage as an mcol
coverageToGRanges <- function(coverage) {
  ends <- cumsum(S4Vectors::runLength(coverage))
  valid.chrs <- lapply(ends, length)>0 #remove the chromosomes with no data
  ends <- ends[valid.chrs] 
  coverage.lvl <- S4Vectors::runValue(coverage)[valid.chrs]
  
  starts <- lapply(ends, function(x) {return(c(1, (x[-length(x)]+1)))})
  
  #convert into a GRanges
  coverage.gr <- lapply(names(ends), 
                        function(chr) {
                          return(toGRanges(
                            data.frame(chr, starts[[chr]], ends[[chr]],
                                       coverage=coverage.lvl[[chr]])
                          )
                          )})
  return(unlist(GenomicRanges::GRangesList(coverage.gr)))
}

#Read a file into a GRanges
fileToGRanges <- function(A, ..., genome=NULL, sep=NULL, comment.char="#") {
  #BED 
  if(IRanges::tolower(tools::file_ext(A))=="bed") {
    #If it's a bed file, try to load it using rtracklayer::import
    gr <- tryCatch(
      expr = {rtracklayer::import(con=A, format = "BED")},
      error = function(e) {
        #If it has failed, it might be because it's a malformed bed file with a
        #header. Raise an informtative error
        
        #read the first line of the file
        ll <- readLines(A, n=1)
        if(hasHeader(ll, sep="\t")) {
          stop("Error in toGRanges when trying to read the BED file \"", A, "\". BED files do can NOT have headers.\n", e)
        } else {
          stop("Error in toGRanges when trying to read the BED file \"", A, "\". Is it a valid BED file?\n", e)
        }
      })
    return(gr)
  }
  
  # #PLINK files
  # NO special treatment is needed. They will be read as data.frame and the 
  # data.frame processing code will identify them.
  #
  # if(IRanges::tolower(tools::file_ext(A)) %in% 
  #    c(".assoc", ".assoc.fisher", ".assoc.dosage",
  #      ".assoc.linear", ".assoc.logistic")) {
  #   
  # } 
  #
  
  
  #GFF
  #TODO
  
  #VCF
  #TODO
  
  #Other files with custom code?
  #TODO
  
  #If we are here, we have a generic text file. 
  num.skip <- firstNonCommentLine(A, comment.char = comment.char) - 1
  
  #check if there are availabe lines in addition to any comment line
  skip.plus.one <- readLines(A, n=num.skip+1)
  skip.plus.one <- skip.plus.one[skip.plus.one!=""] #Remove any empty lines
  if(length(skip.plus.one)==num.skip) return(GRanges())
  
  #Detect header and separator (if needed)
  ll <- readLines(A, n = num.skip + 5)
  ll <- ll[(num.skip+1):length(ll)]
  if(is.null(sep)) sep <- getSeparator(ll)
  has.head <- hasHeader(ll, sep=sep)
  
  #Read the file into a data.frame
  file.cont <- tryCatch(
    expr = {utils::read.table(file = A, header = has.head, sep = sep, skip = num.skip, comment.char = comment.char, ...)},
    error = function(e) {stop("Error in toGRanges when trying to read the file \"", A, "\": ", e)}
  )
  
  #and convert the data.frame into a GRanges using toGRanges itself
  return(tryCatch(
    expr = {toGRanges(file.cont)},
    error = function(e) {stop("Error in toGRanges when building a GRanges with the content of file \"", A, "\": ", e)}
  ))
  
}
