#' Plot Regions
#' 
#' @description
#' Plots sets of regions
#' 
#' @usage plotRegions(x, chromosome, start=NULL, end=NULL, regions.labels=NULL, regions.colors=NULL, ...)
#' 
#' @param x list of objects to be ploted.
#' @param chromosome character or numeric value indicating which chromosome you want to plot.
#' @param start numeric value indicating from which position you want to plot.
#' @param end numeric value indicating to which position you want to plot.
#' @param regions.labels vector indicating the labels for the y axes. It must have the same length as x.
#' @param regions.colors character vector indicating the colors for the plotted regions. It must have the same length as x.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).
#' 
#' @return A plot is created on the current graphics device.
#' 
#' @examples
#' A <- data.frame(chr=1, start=c(1,15,24,40,50), end=c(10,20,30,45,55))
#' 
#' B <- data.frame(chr=1, start=c(2,12,28,35), end=c(5,25,33,43))
#' 
#' plotRegions(list(A,B), chromosome=1, regions.labels=c("A","B"), regions.colors=3:2)
#' 
#' 
#' @export plotRegions
#' 


plotRegions <- function(x, chromosome, start=NULL, end=NULL, regions.labels=NULL, regions.colors=NULL, ...) {
  
  old.scipen <- options()$scipen
  
  options(scipen=999)
  if(!hasArg(chromosome)) stop("chromosome is missing")
  if(!is.null(start) & !is.numeric(start)) stop("start must be numeric")
  if(!is.null(end) & !is.numeric(end)) stop("end must be numeric")
  if(!is.null(regions.colors) & length(regions.colors)!=length(x)) stop("regions.colors must have the same length as x")
  if(!is.null(regions.labels) & length(regions.labels)!=length(x)) stop("regions.labels must have the same length as x")
  if(is.null(regions.colors)) regions.colors<-rep("black",length(x))
  
  
  if(!hasArg(xlab)) xlab <- "position"
  if(!hasArg(ylab)) ylab <- ""
  if(!hasArg(main)) main <- "regions"
  
  #transforming to GRanges and calculating the axes limits for the plot
  maxx <- minx <- rep(0,length(x))
  for(i in 1:length(x)) {
    x[[i]] <- toGRanges(x[[i]])
    x[[i]] <- x[[i]][seqnames(x[[i]])==chromosome,]
    if(!is.null(start)) x[[i]] <- x[[i]][start(x[[i]])>=start,]
    if(!is.null(end)) x[[i]] <- x[[i]][end(x[[i]])<=end,]
    maxx[i] <- max(start(x[[i]]), end(x[[i]]), na.rm=TRUE)
    minx[i] <- min(start(x[[i]]), end(x[[i]]), na.rm=TRUE)
  }
  maxend <- max(maxx)
  minstart <- min(minx)
  
  #plot definition
  plot.new()
  plot.window(xlim=c(minstart, maxend+((maxend-minstart)/10)), ylim=c(0.5, length(x)+0.5), xlab=xlab, ylab=ylab, main=main)
  axis(1, at=round(seq(minstart,maxend,(maxend-minstart)/10),0), labels=round(seq(minstart,maxend,(maxend-minstart)/10),0))
  axis(2, at=c(length(x):1), las=1, labels=regions.labels)
  
  for(reg in 1:length(x)){
    
    reads <- toDataframe(x[[reg]])
    # sort the reads by their start positions
    reads <- reads[order(reads$start),];
    
    # initialise yread: a list to keep track of used y levels
    yread <- c(minstart - 1);
    ypos <- c(); #holds the y position of the ith segment
    
    
    # iterate over segments
    for (r in 1:nrow(reads)){
      read <- reads[r,];
      start <- read$start;
      placed <- FALSE;
      
      # iterate through yread to find the next availible
      # y pos at this x pos (start)
      y <- 1;
      while(!placed){
        
        if(yread[y] < start){
          ypos[r] <- y;
          yread[y] <- read$end + 0.999999;
          placed <- TRUE;
        } 
        
        # current y pos is used by another segment, increment
        y <- y + 1;
        # initialize another y pos if we're at the end of the list
        if(y > length(yread)){
          yread[y] <- minstart-1;
        }
      }
    }
    
    
    # find the maximum y pos that is used to size up the plot
    lengthy <- length(yread)
    ypos <- abs(ypos-(length(unique(ypos))+1))
    reads$ypos <- (ypos + 1)/lengthy
    aux <- length(x)-reg
    reads$ypos <- reads$ypos+aux
    segments(reads$start, reads$ypos, reads$end + 0.999999, reads$ypos, col=regions.colors[[reg]], lwd=3) 
  }
  box(lwd=1.2)
  options(scipen=old.scipen)
  
}

