
#' Summary of permTestResults objects
#' 
#' @method summary permTestResults
#' @return the summary is printed
#' @keywords internal
#' @export
# summary.permTestResults



summary.permTestResults <- function(object, ...) {
  
  
  if(class(object)!="permTestResults")  stop("object must be a permTestResults object")
      
  cat(paste0("Number of permutations: ", object$ntimes, "\n"))
  if(object$ntimes<20) cat("Note: less than 20 permutations might produce unreliable results\n")
  cat("\n")
  cat(paste0("Alternative: ", object$alternative, "\n\n"))
  cat(paste0("Evaluation of the original region set: ", object$observed, "\n\n"))
  cat(paste0("Summary of the evaluation of the permuted region set: \n"))
  print(summary(object$permuted))
  cat("\n\n")
  cat(paste0("Z-score: ", object$zscore, "\n\n"))
  if(object$pval < 0.001 & object$pval >= 0) code <- "***"
  if(object$pval < 0.01 & object$pval >= 0.001) code <- "**"
  if(object$pval < 0.05 & object$pval >= 0.01) code <- "*"
  if(object$pval < 0.1 & object$pval >= 0.05) code <- "."
  if(object$pval <= 1 & object$pval >= 0.1) code <- " "
  cat(paste0("P-value: ", object$pval, " ", code, "\n"))
  cat(paste0("--- \n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"))
    
}
