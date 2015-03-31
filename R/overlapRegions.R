#' Overlap Regions
#' 
#' @description 
#' return overlap between 2 regios set A and B 
#' 
#' @note 
#' The implementation uses when possible the \code{\link{countOverlaps}} function from \code{IRanges} package.
#'  
#' @usage
#' overlapRegions(A, B, colA=NULL, colB=NULL, type="any", min.bases=1, min.pctA=NULL, min.pctB=NULL, get.pctA=FALSE, get.pctB=FALSE, get.bases=FALSE, only.boolean=FALSE, only.count=FALSE, ...)
#' 
#' @param A a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param type
#'\itemize{
#'\item{\code{AinB}: the region in A is contained in a region in B}
#'\item{\code{BinA}: the region in B is contained in A}
#'\item{\code{within}: the region in A or B is contained in a region in the other region set}
#'\item{\code{equal}: the region in A has the same chromosome, start and end as a region in B}
#'\item{\code{AleftB}: the end of the region from A overlaps the beginning of a region in B}
#'\item{\code{ArightB}: the start of a region from A overlaps the end of a region in B}
#'\item{\code{any}: any kind of overlap is returned}
#'}
#' @param colA numeric vector indicating which columns of A the results will contain (default NULL)
#' @param colB numeric vector indicating which columns of B the results will contain (default NULL)
#' @param min.bases numeric minimun number of bp accepted to define a overlap (default 1)
#' @param min.pctA numeric minimun percentage of bases of A accepted to define a overlap (default NULL)
#' @param min.pctB numeric minimun percentage of bases of B accepted to define a overlap (default NULL)
#' @param get.pctA boolean if TRUE add a column in the results indicating the number percentage of A are involved in the overlap (default FALSE) 
#' @param get.pctB boolean if TRUE add a column in the results indicating the number percentage of B are involved in the overlap (default FALSE)
#' @param get.bases boolean if TRUE add in the results the number of overlapped bases (default FALSE)
#' @param only.boolean boolean if TRUE devolve as result a boolean vector containing the overlap state of each regions of A (default FALSE)
#' @param only.count boolean if TRUE devolve as result the number of regions of A overlapping with B
#' @param ... any additional parameter (are there any left?)
#' 
#' @return
#' the default results is a \code{\link{data.frame}} with at least 5 columns "chr" indicating the chromosome of the appartenence of each overlap, "startA", "endA", "startB", "endB", indicating the coordinates of the region A and B for each overlap
#' "type" that describe the nature of the overlap (see arguments "type") eventually other columns can be added (see see arguments "colA", "colB", "get.pctA", "get.pctB", "get.bases") 
#' 
#' 
#' @seealso \code{\link{plotRegions}}, \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{subtractRegions}}, \code{\link{splitRegions}}, \code{\link{extendRegions}}, \code{\link{commonRegions}}, \code{\link{mergeRegions}}, \code{\link{joinRegions}}
#' 
#' @examples
#' A <- data.frame("chr1", c(1, 5, 20, 30), c(8, 13, 28, 40), x=c(1,2,3,4), y=c("a", "b", "c", "d"))
#' 
#' B <- data.frame("chr1", 25, 35)
#' 
#' overlapRegions(A, B)
#' 
#' @export overlapRegions



overlapRegions<-function(A, B, colA=NULL, colB=NULL, type="any", min.bases=1, min.pctA=NULL, min.pctB=NULL, get.pctA=FALSE, get.pctB=FALSE, get.bases=FALSE, only.boolean=FALSE, only.count=FALSE, ...) {
  
  type<-match.arg(type,c("AinB","BinA","within","equal","AleftB","ArightB","any"))
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(B)) stop("B is missing")
  if(!is.numeric(colA) & !is.character(colA) & !is.null(colA)) stop("colA must be a numeric vector")
  if(!is.numeric(colB) & !is.character(colB) & !is.null(colB)) stop("colB must be a numeric vector")
  if(!is.numeric(min.bases)) stop("min.bases must be numeric")
  if(!is.null(min.pctA)) {
    if(!is.numeric(min.pctA) | length(min.pctA)>1) stop("min.pctA must be a numeric between 0 and 100")
    if(min.pctA<0 | min.pctA>100) stop("min.pctA must be a numeric between 0 and 100")
  }
  if(!is.null(min.pctB)) {
    if(!is.numeric(min.pctB)) stop("min.pctB must be a numeric between 0 and 100")
    if(min.pctB<0 | min.pctB>100) stop("min.pctB must be a numeric between 0 and 100")
  }
  if(!is.logical(get.pctA)) stop("get.pctA must be logical")
  if(!is.logical(get.pctB)) stop("get.pctB must be logical")
  if(!is.logical(get.bases)) stop("get.bases must be logical")
  if(!is.logical(only.boolean)) stop("only.boolean must be logical")
  if(!is.logical(only.count)) stop("only.count must be logical")
  if(only.boolean==TRUE & only.count==TRUE) warning("only.count and only.boolean were set to TRUE at the same time. Only.count will be used and a numeric value will be returned.")
  
  
  
  #remember if A or B are dataframes, since it will affect the way the additional columns are used
  A.was.dataframe = (is.data.frame(A))
  B.was.dataframe = (is.data.frame(B))
  
  A <- toGRanges(A)
  B <- toGRanges(B)
  
  
  ##### If it is possible to directly use a GRanges function and return early, do it #####
  
  if((only.count==TRUE || only.boolean==TRUE) && is.null(min.pctA) && is.null(min.pctB)) { #No limiting by the percentage of overlap (no supported by GRanges)
                                                                                          # and only interested in counts or in boolean answer
    granges.supported.types <- c("any", "within", "equal")
    if(type %in% granges.supported.types) { #using one of the overalp types supported by GRanges
      #Then, use GRanges "directly"  
      ov <- countOverlaps(A, B, type=type, minoverlap=min.bases)
      if(only.count) {
        return(sum(ov))
      }
      if(only.boolean) {
        return(ov>0)
      }
    }
  }
  
  
  ####  If it's not possible to directly use GRanges filters, do it ourselves   ####
  
  #Find the overlaps (partially filtered)
    #use the findOverlaps filtering as much as possible
    ty <- type
    if(ty=="AinB") {
      ty <- "within"
    }
    if(ty=="AleftB" || ty=="ArightB" || ty=="BinA") {
      ty <- "any"
    }
  
    over.AB<-findOverlaps(A, B, type=ty, minoverlap=min.bases)
  
  
    tab.A<-A[queryHits(over.AB)]
    tab.B<-B[subjectHits(over.AB)]
  
  
    #Add the additional columns specified in colA and colB
      #if A or B were dataframes and colA or colB are numeric, we have to substract 3 to get the right columns,
      #since the first three columns are not included in mcols when transformed into GRanges
      if(!is.null(colA) && A.was.dataframe && is.numeric(colA)) {colA <- colA -3}
      if(!is.null(colB) && B.was.dataframe && is.numeric(colB)) {colB <- colB -3}
  
      #now colA corresponds to the columns in mcols
  
    
  
  #Build the results table. If its empty, build it anyway
  if(length(tab.A)>0) {
    tab<-data.frame(chr=as.character(seqnames(tab.A)), startA=start(tab.A), endA=end(tab.A), startB=start(tab.B), endB=end(tab.B))
  } else {
    tab<-data.frame(chr=character(), startA=numeric(), endA=numeric(), startB=numeric(), endB=numeric())
  }

  #names(tab) <- c("chr", "startA", "endA", "startB", "endB") #Somehow the name of the first column was changed to "value", set the names again

  
  if(!is.null(colA)) {
    if(length(colA)==1 && colA=="all") {
      dd <- as.data.frame(mcols(tab.A))
      col.with.name <- 1
    } else {
      dd <- as.data.frame(mcols(tab.A))[,colA]
      col.with.name <- colA
    }
    tab <- cbind(tab, dd)
    if(length(mcols(tab.A))==1) { #If the length is only one, the column name is "lost in conversion"
      names(tab)[length(tab)] <- names(mcols(tab.A))[col.with.name]  #So, recover it
    }
  }
  
  if(!is.null(colB)) {
    if(length(colB)==1 && colB=="all") {
      dd <- as.data.frame(mcols(tab.B))
      col.with.name <- 1
    } else {
      dd <- as.data.frame(mcols(tab.B))[,colB]
      col.with.name <- colB
    }
    tab <- cbind(tab, dd)
    if(length(mcols(tab.B))==1) { #If the length is only one, the column name is "lost in conversion"
      names(tab)[length(tab)] <- names(mcols(tab.B))[col.with.name]  #So, recover it
    }
  }
  
  #Ordering is not needed, since the result from findOvelaps is already ordered
    #tab<-tab[order(tab[,1]),]
  
  #Compute the overlap type for all the overlaps
  vec.type<-vector(mode="character", length=dim(tab)[1])
  vec.type[which((tab$startA >= tab$startB) & (tab$endA <= tab$endB))]<-"AinB"
  vec.type[which((tab$startA <= tab$startB) & (tab$endA >= tab$endB))]<-"BinA"
  vec.type[which((tab$startA == tab$startB) & (tab$endA == tab$endB))]<-"equal"
  vec.type[which((tab$startA < tab$startB) & (tab$endA < tab$endB))]<-"AleftB"
  vec.type[which((tab$startA > tab$startB) & (tab$endA > tab$endB))]<-"ArightB"
  
  tab <- cbind(tab, type=vec.type)
 
  
  #Filter out any remaining unwanted overlap  
  if(type=="AinB") tab<-tab[tab[,"type"]=="AinB" | tab[,"type"]=="equal",]
  if(type=="BinA") tab<-tab[tab[,"type"]=="BinA" | tab[,"type"]=="equal",]
  if(type=="AleftB") tab<-tab[tab[,"type"]=="AleftB",]
  if(type=="ArightB") tab<-tab[tab[,"type"]=="ArightB",]
  

                              
  #Apply overlap size filters and add this information if requested
  if(any(!is.null(min.pctA), !is.null(min.pctB), get.pctA, get.pctB, get.bases)) {                             
                        
  
    start<-apply(cbind(tab$startA, tab$startB), 1, max)
    end<-apply(cbind(tab$endA, tab$endB), 1, min)
    ov.bases <- end-start +1
      
    if(get.bases==TRUE) tab <- cbind(tab, ov.bases)
 
    pct.basesA <- (ov.bases/(tab$endA-tab$startA))*100.0
    pct.basesB <- (ov.bases/(tab$endB-tab$startB))*100.0
    
    if(get.pctA==TRUE) tab <- cbind(tab, pct.basesA)
    if(get.pctB==TRUE) tab <- cbind(tab, pct.basesB)
    
    
    to.keep <- as.logical(rep(TRUE, length(tab[,1])))
    if(!is.null(min.pctA)) to.keep <- to.keep & (pct.basesA >= min.pctA)
    if(!is.null(min.pctB)) to.keep <- to.keep & (pct.basesB >= min.pctB)
    
    tab <- tab[to.keep, ]
    
  }
  
  #Certain combinations of parameters (e.g. asking for min.pct filtering) will lead us here even if only.boolean or only.count were set.
  #So take care of these cases.
  if(only.count) {
    return(length(tab[,1])) 
  }
  if(only.boolean) {
    #for every region in A: check if it's in tab
    #To do it, build an id for each region in A and tab using chr, start and end and compare the ids
    #Could be done with a call to overlap with type "equal" but I think it would be slower. TODO: Check
    A.ids <- paste(seqnames(A), start(A), end(A), sep="#")
    tab.ids <- paste(tab[,1], tab[,2], tab[,3], sep="#")
    return(A.ids %in% tab.ids)    
  }
  
  return(tab)
}

