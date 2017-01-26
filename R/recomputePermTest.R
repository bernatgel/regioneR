#' Recompute Permutation Test
#' 
#' @description
#' Recomputes the permutation test changing the alternative hypotesis
#' 
#' @usage recomputePermTest(ptr)
#' 
#' @param ptr   an object of class \code{permTestResults}
#'      
#' @return
#' A list of class \code{permTestResults} containing the same components as \code{\link{permTest}} results.
#' 
#' @seealso \code{\link{permTest}}
#' 
#' @examples
#' A <- createRandomRegions(nregions=10, length.mean=1000000)
#' 
#' B <- createRandomRegions(nregions=10, length.mean=1000000)
#' 
#' resPerm <- permTest(A=A, B=B, ntimes=5, alternative="less", genome="hg19", evaluate.function=meanDistance, randomize.function=randomizeRegions)
#' 
#' plot(resPerm)
#' 
# resPermRecomputed <- recomputePermTest(resPerm)
#' 
# summary(resPermRecomputed)
#' 
# plot(resPermRecomputed)
#'  
#' @export recomputePermTest



recomputePermTest<-function(ptr){
  
  if(class(ptr)!="permTestResults")  stop("x must be a permTestResults object")
  
  ptr2<-ptr
  if(ptr$alternative == "less"){
    ptr2$pval <- (sum(ptr$observed <= ptr$permuted) + 1) / (ptr$ntimes + 1)
    ptr2$alternative<-"greater"
  } 
  if(ptr$alternative == "greater"){
    ptr2$pval <- (sum(ptr$observed >= ptr$permuted) + 1) / (ptr$ntimes + 1)
    ptr2$alternative<-"less"
  } 
  return(ptr2) 
}
