#' Split Regions
#' 
#' @description
#' Splits a region set A by both ends of the regions in a second region set B.
#' 
#' @usage splitRegions(A, B, min.size=1, track.original=TRUE)
#' 
#' @param A a region set in any of the formats accepted by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B a region set in any of the formats accepted by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param min.size numeric value, minimal size of the new regions
#' @param track.original logical indicating if you want to keep the original regions and additional information in the output
#' 
#' @return
#' A GRanges with the splitted regions.
#' 
#' @seealso \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{subtractRegions}}, \code{\link{commonRegions}}, \code{\link{extendRegions}}, \code{\link{joinRegions}}, \code{\link{mergeRegions}}, \code{\link{overlapRegions}}
#' 
#' @examples
#' A <- data.frame(chr=1, start=c(1, 15, 24, 40, 50), end=c(10, 20, 30, 45, 55))
#' 
#' B <- data.frame(chr=1, start=c(2, 12, 28, 35), end=c(5, 25, 33, 43))
#' 
#' splits <- splitRegions(A, B)
#' 
#' plotRegions(list(A, B, splits), chromosome=1, regions.labels=c("A", "B", "splits"), regions.colors=3:1)
#' 
#' @export splitRegions


splitRegions <- function(A, B, min.size=1, track.original=TRUE) {

  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")
  if(!is.numeric(min.size)) stop("min.size must be numeric")
  if(!is.logical(track.original)) stop("track.original must be logical")
  
  
  A <- toDataframe(A)
  B <- toDataframe(B)
  
#creating the outputfiles
  if(track.original == FALSE) outA <- rep(NA, 3)
  if(track.original == TRUE) outA <- rep(NA, dim(A)[2]+2)
  
#for each chromosome
  chr<-unique(c(as.character(as.vector(A[,1])), as.character(as.vector(B[,1]))))
  for(i in 1:length(chr)) {
    Achr <- A[A[,1]==chr[i],]
    Bchr <- B[B[,1]==chr[i],]
    

	#adding all the starts and ends of A and B and indicating if they are start (s) or end (e)
    splits <- data.frame(x=c(Achr[,2], Bchr[,2], Achr[,3], Bchr[,3]),y=c(rep("s", dim(Achr)[1]), 
                        rep("s", dim(Bchr)[1]),rep("e", dim(Achr)[1]),rep("e", dim(Bchr)[1])), 
                        stringsAsFactors=FALSE)
    #sorting
    splits <- splits[order(splits[,1]),]
    #removing duplicated rows
    #splits<-splits[!duplicated(splits),]
    
    #once we have the splits we want to find the regions
    start <- end <- 0
    for(j in 1:(dim(splits)[1]-1)) {
      #if the split is start, the start point of the region will be the split but the end point of the region will be th split minus 1
      #if the split is end,the start point of the region will be the split plus 1 and the end point of the region will be the next split
      if(splits[j,2]=="s") {
        start <- c(start,splits[j,1])
      } else {
        start <- c(start,splits[j,1]+1)
      }
      if(splits[(j+1),2]=="s") {
        end <- c(end,splits[(j+1),1]-1)
      } else {
        end <- c(end,splits[j+1,1])
      }
      
    }
    start <- start[-1]
    end <- end[-1]
    newregions <- data.frame(start,end)
    #removing regions with no length
    long <- end-start+1
    newregions <- newregions[long>0,]
    
    #constructing the output files and adding the original info if track.original==TRUE
    for(j in 1:dim(newregions)[1]) {
      auxA <- which(Achr[,2]<=newregions[j,1] & Achr[,3]>=newregions[j,2])
      if(length(auxA)>0) {
        if(track.original==TRUE){
          d1 <- newregions[j,]
          d2 <- d1[rep(seq_len(nrow(d1)), length(auxA)), ]
          outA <- rbind(outA, data.frame(Achr[auxA,1], d2, Achr[auxA,2:dim(Achr)[2]]))
        }
          if(track.original==FALSE) outA <- rbind(outA, c(chr[i], newregions[j,]))
      }
    }
  }

  colnames(outA)[1] <- "chr"
  rownames(outA) <- NULL
  outA <- outA[-1,]
    
  res <- toGRanges(outA)
  return(res)
}





