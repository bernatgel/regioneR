
#' Summary of permTestResultsList objects
#' 
#' @method summary permTestResultsList
#' @return the summary is printed
#' @keywords internal
#' @export
# summary.permTestResultsList



summary.permTestResultsList <- function(object, ...) {
  
  
  if(class(object)!="permTestResultsList")  stop("object must be a permTestResultsList object")
  
  if(length(object) > 0) {
    
    res <- do.call(rbind, lapply(object, function(x) {
      return(data.frame(pvalue=x$pval, zscore=x$zscore, test=x$alternative))
    }))
    
    
    cat(paste0("Permutation tests: ", length(object), "\n"))
    cat(paste0("Significant permutation tests: ", length(which(res$pvalue<=0.05)), "\n"))
    cat(paste0("Iterations: ", object[[1]]$ntimes, "\n"))
    cat(paste0("Randomization Function: ", object[[1]]$randomize.function.name, "\n"))
    cat("Tests Results:\n")
    print(res)
  } else {
    cat(paste0("permTestResultList object of length 0\n"))
  }
     
  
}

