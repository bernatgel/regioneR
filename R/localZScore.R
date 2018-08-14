#' Local z-score
#' 
#' @description 
#' Evaluates tthe variation of the z-score in the vicinty of the original region set
#' 
#' @usage 
#' localZScore(A, pt, window, step, ...)
#'  
#' @param A a region set in any of the formats accepted by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param pt a permTestResult object
#' @param window a window in wich the local Z-score will be calculated (bp)
#' @param step the number of bp that divide each Z-score evaluation  
#' @param ... further arguments to be passed to other methods.
#' 
#' @return
#' It returns a local z-score object
#' 
#' @seealso \code{\link{overlapPermTest}}, \code{\link{permTest}}
#' 
#' @examples
#' 
#' genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
#' A <- createRandomRegions(nregions=20, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
#' B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
#' 
#' pt <- overlapPermTest(A=A, B=B, ntimes=10, genome=genome, non.overlapping=FALSE)
#' plot(pt)
#'  
#' lz <- localZScore(A=A, B=B, pt=pt)
#' plot(lz)
#' 
#' 
#' pt2 <- permTest(A=A, B=B, ntimes=10, randomize.function=randomizeRegions, evaluate.function=list(overlap=numOverlaps, distance=meanDistance), genome=genome, non.overlapping=FALSE)
#' plot(pt2)
#' 
#' lz2 <- localZScore(A=A, B=B, pt2)
#' plot(lz2)
#' 
#'  
#' @export localZScore
#' 
#' @importFrom GenomicRanges shift
  





localZScore <- function(A, pt, window, step, ...) {
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(pt)) stop("pt is missing")
  
  A <- toGRanges(A)
  
  if(!hasArg(window)) {
    window <- 5*mean(width(A))
  }
  if(!hasArg(step)) {
    step <- floor(window/10)
  }
  
  #if pt is a permTestResults object, compute the localZCScore and return it
  if(class(pt) == "permTestResults") {
    mean.permuted <- mean(pt$permuted)
    sd.permuted <- stats::sd(pt$permuted)
    
    num.steps <- floor(window/step)
    
    shifts <- (1:num.steps)*step
    shifts <- c(rev(-1*shifts), 0, shifts)
    
    shifted.z.score <- function(shift) {
      shifted.A <- GenomicRanges::shift(A, shift)
      shifted.evaluation <- tryCatch(pt$evaluate.function(shifted.A, ...),
                                     error=function(e) {stop(paste0("There was an error when computing evaluation function of the shifted region set: \n", as.character(e),
                                                             "Evaluation Function: ", pt$evaluate.function.name,
                                                             "\n Shift: ", shift))}
                                      )
      shifted.z.score <- (shifted.evaluation - mean.permuted)/sd.permuted
      return(shifted.z.score)
    }
    
    shifted <- lapply(as.list(shifts), shifted.z.score)
    
    shifted.z.scores <- do.call(c, shifted)
    
    rZ <- list(shifted.z.scores=shifted.z.scores, shifts=shifts, window=window, step=step, original.z.score=pt$zscore)
    class(rZ) <- "localZScoreResults"
    return(rZ)
    
  } else { #else, if it's a list, run localZscore for each element in the list
    if(class(pt) == "permTestResultsList") {
      lz <- mclapply(pt, function(ptt) {return(localZScore(A=A, pt=ptt, window=window, step=step, ...))})
      class(lz) <- "localZScoreResultsList"
      return(lz)
    } else {
      stop(paste0("pt must be of class permTestResults or permTestResultsList and is a: ", class(pt)))
    }
  }
}

