#' Overlap Graphical Summary
#' 
#' @description
#' Graphical summary of the overlap between two set of regions.
#' 
#' @usage overlapGraphicalSummary(A, B, regions.labels=c("A","B"), regions.colors=c("black","forestgreen","darkred"), ...)
#' 
#' @param A a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param B a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param regions.labels vector indicating the labels for the y axes.
#' @param regions.colors character vector indicating the colors for the regions.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}). 
#'
#'  @return A plot is created on the current graphics device.
#'   
#' @seealso \code{\link{overlapPermTest}}, \code{\link{overlapRegions}}
#' 
#' @examples
#' A <- data.frame(chr=1, start=c(1,15,24,40,50), end=c(10,20,30,45,55))
#' 
#' B <- data.frame(chr=1, start=c(2,12,28,35), end=c(5,25,33,43))
#' 
#' overlapGraphicalSummary(A, B, regions.labels=c("A","B"), regions.colors=c(4,5,6))
#' 
#' @export overlapGraphicalSummary


overlapGraphicalSummary <- function(A, B, regions.labels=c("A","B"), regions.colors=c("black","forestgreen","darkred"), ...) {
  
  if(!hasArg(A)) stop("A is required")
  if(!hasArg(B)) stop("B is required")
  
  if(!hasArg(main)) main <- ""
  
  old.scipen <- options()$scipen
  
  options(scipen=999)
  
  A <- toGRanges(A)
  B <- toGRanges(B)
  
  test <- overlapRegions(A, B, get.pctA=TRUE, get.pctB=TRUE, get.bases=TRUE)
  test.gr <- toGRanges(test)
  
  lin <- (test$pct.basesA + test$pct.basesB) / 2
  linA <- (test$pct.basesA)
  linB <- (test$pct.basesB)
  
  layout(matrix(c(1,2,3,4,5,6), nrow=2, byrow = TRUE), widths=c(3,1.5,3), heights=c(3,3))
  pl <- plot(lin[order(lin)], type="l", lwd=5, ylim=c(0,110), main="", col=regions.colors[1], las=1, ylab="", xlab="", ...)
  abline(h=50, lty=2)
  abline(h=100, lty=2)
  lines(linA[order(linA)], col=regions.colors[2], lwd=3)
  lines(linB[order(linB)], col=regions.colors[3], lwd=3)
  title(paste("overlap ",regions.labels[1]," on ",regions.labels[2]," = ",round((length(unique(test.gr))/length(A))*100, digits=2),"%"))
  box(lwd=1.2)
  
  par(xpd=TRUE)
  plot(1, type="n", axes=FALSE, xlab="", ylab="", ...)
  legend("right", c(paste0(regions.labels[1],"+",regions.labels[2]),paste0(regions.labels)), col=regions.colors, lwd=3)
  par(xpd=FALSE)
  
  C <- c(length(A), length(B), length(unique(test.gr)), sum(test$pct.basesA==100), sum(test$pct.basesB==100), sum((test$pct.basesA==100)&(test$pct.basesB==100)))
  namearg <- c(paste0("reg.",regions.labels[1]),paste0("reg.",regions.labels[2]),paste0("uni.ov.",regions.labels[1]),paste0(regions.labels[1],".in.",regions.labels[2]),paste0(regions.labels[2],".in.",regions.labels[1]),paste0(regions.labels[1],"=",regions.labels[1]))
  barplot(C, main="overlap types", names.arg=namearg, las=2, ylim=c(0, max(C)+0.1), ...)
  box(lwd=1.2)
  
  pieA <- c(length(A)-length(unique(test.gr)), length(A)-(length(A)-length(unique(test.gr))))
  pie(pieA, labels = paste0(pieA), main = "Overlap of A", col=c("white","gray"))
  box(lwd=1.2)
  
  par(xpd=TRUE)
  plot(1, type="n", axes=FALSE, xlab="", ylab="", ...)
  legend("center", legend=c("non overlaped","overlaped"), fill=c("white","gray"))
  par(xpd=FALSE)
  
  pieB<-c(length(B)-length(unique(test.gr)),length(B)-(length(B)-length(unique(test.gr))))
  pie(pieB, labels = paste0(pieB), main = "Overlap of B", col=c("white","gray"))
  box(lwd=1.2)
  
  
  par(mfrow=c(1,1))
  options(scipen=old.scipen)
  
}




