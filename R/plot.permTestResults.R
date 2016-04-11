# Plot Permutation Test Results
#
# @description
#' Function for plotting the results from a \code{permTestResults} object.
#' 
#' @method plot permTestResults
#' 
#' @param x          an object of class \code{permTestResults}.
#' @param pvalthres  p-value threshold for significance. Default is 0.05.
#' @param plotType   the type of plot to display. This must be one of \code{"Area"} or \code{"Tailed"}. Default is \code{"Area"}.
#' @param main       a character specifying the title of the plot. Defaults to "".
#' @param xlab       a character specifying the label of the x axis. Defaults to NULL, which produces a plot with the evaluation function name as the x axis label.
#' @param ylab       a character specifying the label of the y axis. Defaults to "".
#' @param ...        further arguments to be passed to or from methods.
#' 
#' @return A plot is created on the current graphics device.
#' 
#' @seealso \code{\link{permTest}}
#' 
#' @examples 
#' 
#' genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
#' A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
#' B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
#' 
#' pt <- overlapPermTest(A=A, B=B, ntimes=10, genome=genome, non.overlapping=FALSE)
#' summary(pt)
#' plot(pt)
#' plot(pt, plotType="Tailed")  
#'  
#' pt2 <- permTest(A=A, B=B, ntimes=10, alternative="auto", genome=genome, evaluate.function=meanDistance, randomize.function=randomizeRegions, non.overlapping=FALSE)
#' summary(pt2)
#' plot(pt2)
#' plot(pt2, plotType="Tailed")
#' 
#' @import graphics
#' @importFrom stats dnorm qnorm rnorm runif
#'
#' @export 


plot.permTestResults<-function(x, pvalthres=0.05, plotType="Tailed", main="", xlab=NULL, ylab="", ...){
  
  old.scipen <- options()$scipen
  
  options(scipen=999)
  
  if(class(x)!="permTestResults")  stop("x must be a permTestResults object")
  if(!is.numeric(pvalthres)) stop("pvalthres must be numeric")
  plotType<-match.arg(plotType,c("Area","Tailed"))
  
  
  if(is.null(xlab)) xlab <- paste0(x$evaluate.function.name)
  if(nchar(main)>0) main <- paste0(main, "\n")  
  
  alternative<-x$alternative
  xcoords<-x$permuted
  xcoords<-xcoords[order(xcoords)]
  pval<-round(x$pval,4)
  nperm<-x$ntimes
  mperm<-mean(xcoords,na.rm=TRUE)
  mobs<-x$observed
  zscore<-round(x$zscore,3)
  
  if(is.finite(zscore)){
    y<-dnorm(xcoords,mean=mean(xcoords,na.rm=TRUE),sd=sd(xcoords,na.rm=TRUE))
    xhist<-hist(xcoords,breaks=30,plot=FALSE)$density
    ymax<-max(max(y,na.rm=TRUE),max(xhist,na.rm=TRUE))
    
    if (alternative=="greater") aux<-qnorm((1-pvalthres),mean=mean(xcoords,na.rm=TRUE),sd=sd(xcoords,na.rm=TRUE))
    if (alternative=="less") aux<-qnorm(pvalthres,mean=mean(xcoords,na.rm=TRUE),sd=sd(xcoords,na.rm=TRUE))
    
    xmin<-min(mobs,min(xcoords,na.rm=TRUE),min(aux,na.rm=TRUE),na.rm=TRUE)
    xmax<-max(mobs,max(xcoords,na.rm=TRUE),max(aux,na.rm=TRUE),na.rm=TRUE)
    
    hist(xcoords, prob = TRUE, ylim = c(0,ymax), breaks = 30, xlim = c(xmin,xmax),
         las = 1, col = "lightgray", border = "lightgray", xlab=xlab, ylab=ylab, main=paste(main, "p-value: ",
                                                                                            pval, "\n Z-score: ", zscore, "\n n perm: ", nperm, "\n randomization: ", paste0(x$randomize.function.name)),
         cex.main=0.8, ...)
    
    if(plotType=="Area"){
      if(alternative=="greater"){
        polygon(c(aux,aux,xmax,xmax),c(max(y,na.rm=TRUE),0,0,max(y,na.rm=TRUE)),col="red",density=10,border="white")
        lines(c(aux,aux),c(0,ymax*0.8),col="red",lwd=3)
        text(aux,ymax*0.9,bquote(alpha==.(pvalthres)),cex=0.8,pos=4)
      }
      if(alternative=="less"){
        polygon(c(aux,aux,xmin,xmin),c(max(y,na.rm=TRUE),0,0,max(y,na.rm=TRUE)),col="red",density=10,border="white")
        lines(c(aux,aux),c(0,ymax*0.8),col="red",lwd=3)
        text(aux,ymax*0.9,bquote(alpha==.(pvalthres)),cex=0.8)
      }
    }
    
    if(plotType=="Tailed"){
      if(alternative=="greater"){
        aux3<-seq(aux,xmax,length=50)
        y3<-dnorm(aux3,mean(xcoords,na.rm=TRUE),sd(xcoords,na.rm=TRUE))
        polygon(c(aux3[1],aux3,aux3[length(aux3)]),c(0,y3,0),col="red",density=10,border="white")
        lines(c(aux,aux),c(0,ymax*0.8),col="red",lwd=3)
        text(aux,ymax*0.9,bquote(alpha==.(pvalthres)),cex=0.8)
      }
      if(alternative=="less"){
        aux3<-seq(aux,xmin,length=50)
        y3<-dnorm(aux3,mean(xcoords,na.rm=TRUE),sd(xcoords,na.rm=TRUE))
        polygon(c(aux3[1],aux3,aux3[length(aux3)]),c(0,y3,0),col="red",density=10,border="white")
        lines(c(aux,aux),c(0,ymax*0.8),col="red",lwd=3)
        text(aux,ymax*0.9,bquote(alpha==.(pvalthres)),cex=0.8)
      }
    }
    
    
    lines(xcoords,y,lwd=2)
    lines(c(mperm,mperm),c(0,ymax*0.8),col="black",lwd=3)
    text(mperm,ymax*0.9,expression(Ev[perm]),cex=0.8)
    
    lines(c(mobs,mobs),c(0,ymax*0.8),col="forestgreen",lwd=3)
    text(mobs,ymax*0.9,expression(Ev[obs]),cex=0.8)
    arrows(mperm,ymax*0.75,mobs,ymax*0.75,length=0.1,code=3)
    box(lwd=1.2)
  }
  
  if(!is.finite(zscore)){
    xhist<-hist(xcoords,breaks=30,plot=FALSE)$density
    ymax<-max(xhist,na.rm=TRUE)
    
    hist(xcoords, prob = TRUE, ylim = c(0,ymax), breaks = 30, xlim = c(min(mobs,min(xcoords,na.rm=TRUE),na.rm=TRUE), max(mobs,max(xcoords,na.rm=TRUE),na.rm=TRUE)), las = 1, col = "lightgray", border = "lightgray", xlab=xlab, ylab=ylab, main=paste(main, "p-value: ", pval, "\n Z-score: ", zscore, "\n n perm: ", nperm, "\n randomization: ", paste0(x$randomize.function.name)), cex.main=0.8, ...)
    
    lines(c(mperm,mperm),c(0,ymax*0.8),col="black",lwd=3)
    text(mperm,ymax*0.9,expression(Ev[perm]),cex=0.8)
    
    lines(c(mobs,mobs),c(0,ymax*0.8),col="forestgreen",lwd=3)
    text(mobs,ymax*0.9,expression(Ev[obs]),cex=0.8)
    arrows(mperm,ymax*0.75,mobs,ymax*0.75,length=0.1,code=3)
    box(lwd=1.2)
    
    
    warning(paste0("all permuted values are equal to ",xcoords[1],". It is not posible to adjust a normal distribution nor to compute a Z-score."))
    
  }
  
  options(scipen=old.scipen)
  
}
