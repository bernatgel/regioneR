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
#' @usage toGRanges(A, ...)
#' 
#' @param A  a \code{\link{data.frame}} containing a region set, a \code{\link{GRanges}} object, a BED file or any type of file supported by \code{rtracklayer}
#' @param ... further arguments to be passed to other methods.
#' 
#' @return
#' A \code{\link{GRanges}} object with the regions in A
#' 
#' @seealso \code{\link{toDataframe}}
#' 
#' @examples
#' A <- data.frame(chr=1, start=c(1, 15, 24), end=c(10, 20, 30),  x=c(1,2,3), y=c("a", "b", "c"))
#' 
#' toGRanges(A)
#' 
#' @export toGRanges
#' 
#' @import GenomicRanges
#' @import BSgenome
#' @import memoise
#' @importFrom rtracklayer import
#' @import parallel
#' @import GenomeInfoDb
#' @import IRanges
#' 
#' 

#' 


#STRANDED VERSION
#In addition, the function will check if there's a column 
#named "strand" and use it to set the strand information of the GRanges object. If no "strand" column is present, if a fourth column 
#exists and it contains only "+", "-", and "*", that column will be used as strand.
# 
# toGRanges <- function(A, stranded=TRUE) {
#   
#   if(!hasArg(A)) stop("A is missing")
#   
#   if(is(A, "GRanges")) {
#     return(A)
#   }
#   
#   
#   #if it's a dataframe, assume the three first columns are: chr, start and end
#   if(is(A, , "data.frame")) {
#     gr<-GRanges(seqnames=A[,1],ranges=IRanges(A[,2],end=A[,3]))
#     #We cannot assign the metadata in a single line because when only one metadata column was requested, 
#     #it was automatically transformed into a vector and it lost its name. 
#     if(dim(A)[2]>3) { #if there's metadata or strand information
#       
#       if(stranded) {
#         if(!is.null(A$strand)) {
#           strand.col <- which(names(A)=="strand")[1]
#         } else {
#           strand.col <- grep("strand", names(A))  
#           if(length(strand.col)>0) {
#             strand.col <- strand.col[1]        
#           } else {
#             strand.col <- 4
#           }
#         }    
#         
#         
#         if(all(names(table((A[,strand.col]))) %in% c("+", "-", "*"))) { #this might slower than using "levels", but it will work even if it's not a factor but a character.
#           strand(gr) <- A[,strand.col]
#         } else {
#           strand.col<-NULL
#         }
#       } else {
#         strand.col <- NULL
#         #Change the name of the strand column, if any, so it does not interfere with the GRanges creation process
#         if(!is.null(A$strand)) {
#           names(A)[names(A)=="strand"] <- "strand.1"
#           #QUESTION: Do we want to raise a warning?
#         }
#       }
#       
#       #Store the metadata
#       metadata.columns <- c(4:ncol(A))
#       if(!is.null(strand.col)) { #Remove the strand column if present
#         metadata.columns <- metadata.columns[-c(strand.col-3)]
#       }
#       if(length(metadata.columns>0)) {
#         original.metadata <- as.data.frame(A[,metadata.columns])   
#         names(original.metadata) <- names(A)[metadata.columns]
#         elementMetadata(gr) <- original.metadata
#       }
#       
#     }
#     return(gr)
#   }
#   
#   #If its something else, assume it's something rtracklayer can read and import from it
#   gr <- rtracklayer::import(con=A)
#   return(gr)  
#   
# }


toGRanges <- function(A, ...) {
    
  if(!hasArg(A)) stop("A is missing")
  
  if(is(A, "GRanges")) {
    return(A)
  }
  
  gr <- NULL
  
  #if it's a dataframe, assume the three first columns are: chr, start and end
  if(is(A, "data.frame")) {
    
    chrs <- as.character(A[,1]) #Transform the first column into a character. It does not work with factors.
    
    gr<-GRanges(seqnames=chrs,ranges=IRanges(A[,2],end=A[,3]))
    #We cannot assign the metadata in a single line because when only one metadata column was requested, 
    #it was automatically transformed into a vector and it lost its name. 
    if(dim(A)[2]>3) { #if there's metadata or strand information
      
    
      #Store the metadata
      metadata.columns <- c(4:ncol(A))
      if(length(metadata.columns>0)) {
        original.metadata <- as.data.frame(A[,metadata.columns])   
        names(original.metadata) <- names(A)[metadata.columns]
        elementMetadata(gr) <- original.metadata
      }
      
    }
    return(gr)
  }
  
  #If its something else, try with rtracklayer and see if it can read and import from it.
  #If it fails and raises an error (usually because unknown format (txt)), try to open it with read.delim and if it fails, try with read.csv
  gr <- tryCatch(
    expr = {return(rtracklayer::import(con=A, ...))},
    error = function(err) {
      return(tryCatch(
        expr = {return(toGRanges(read.delim(A, ...)))},
        error = function(err) {
          return(toGRanges(read.csv(A, ...)))
        }
      ))
    }
  )

  return(gr)  
  
}
