#' Print permTestResults objects
#' @return the object is printed
#' 
#' @examples
#' genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
#' A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
#' B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
#' 
#' pt <- permTest(A=A, B=B, ntimes=10, alternative="auto", verbose=TRUE, genome=genome, evaluate.function=meanDistance, randomize.function=randomizeRegions, non.overlapping=FALSE)
#' print(pt)
#' 
#' @keywords internal
#' @export 


print.permTestResults <- function(x, ...) {
  cat(paste0("P-value: ", x$pval, "\n"))
  cat(paste0("Z-score: ", x$zscore, "\n"))
  cat(paste0("Number of iterations: ", x$ntimes, "\n"))
  if(x$ntimes<20) cat("Note: less than 20 iterations might produce unreliable results\n")
  cat(paste0("Alternative: ", x$alternative, "\n"))
  cat(paste0("Evaluation of the original region set: ", x$observed, "\n"))
  cat(paste0("Evaluation function: ", x$evaluate.function.name, "\n"))
  cat(paste0("Randomization function: ", x$randomize.function.name, "\n"))
}



