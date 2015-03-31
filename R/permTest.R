#' Permutation Test
#' 
#' Performs a permutation test to see if there is an association between a region set and some other feature using
#' an evaluation function.
#' #' 
#' @usage permTest(A, ntimes=100, randomize.function, evaluate.function, alternative="auto", min.parallel=1000, force.parallel=NULL, randomize.function.name=NULL, evaluate.function.name=NULL, verbose=FALSE, ...)
#' 
#' @param A a region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param ntimes number of permutations
#' @param randomize.function function to create random regions. It must return a set of regions.
#' @param evaluate.function function to search for association. It must return a numeric value.
#' @param alternative the alternative hypothesis must be one of \code{"greater"}, \code{"less"} or \code{"auto"}. If \code{"auto"}, the alternative will be decided depending on the data.
#' @param min.parallel if force.parallel is not specified, this will be used to determine the threshold for parallel computation. If \code{length(A) * ntimes > min.parallel}, it will activate the parallel computation. Single threaded otherwise.
#' @param force.parallel logical indicating if the computation must be paralelized. 
#' @param randomize.function.name character. If specified, the permTestResults object will have this name instead of the name of the randomization function used. Useful specially when using unnamed anonymous functions.
#' @param evaluate.function.name character. If specified, the permTestResults object will have this name instead of the name of the evaluation function used. Useful specially when using unnamed anonymous functions.
#' @param verbose a boolean. If verbose=TRUE it creates a progress bar to show the computation progress. When combined with parallel computation, it might have an impact in the total computation time.
#' @param ... further arguments to be passed to other methods.
#' 
#' @details permTest performs a permutation test of the regions in RS to test the association with the feature evaluated with the evaluation function.
#' The regions are randomized using the randomization.function and the evaluation.function is used to evaluate them. More information can be found in
#'  the vignette.
#' 
#' @return
#' A list of class \code{permTestResults} containing the following components:
#' \itemize{
#' \item \bold{\code{pval}} the p-value of the test.
#' \item \bold{\code{ntimes}} the number of permutations.
#' \item \bold{\code{alternative}} a character string describing the alternative hypotesis.
#' \item \bold{\code{observed}} the value of the statistic for the original data set.
#' \item \bold{\code{permuted}} the values of the statistic for each permuted data set.
#' \item \bold{\code{zscore}} the value of the standard score. \code{(observed-\link{mean}(permuted))/\link{sd}(permuted)}
#' \item \bold{\code{randomize.function}} the randomization function used.
#' \item \bold{\code{randomize.function.name}} the name of the randomization used.
#' \item \bold{\code{evaluate.function}} the evaluation function used.
#' \item \bold{\code{evaluate.function.name}} the name of the evaluation function used.
#' }
#' 
#' @seealso \code{\link{overlapPermTest}}
#' 
#' @examples 
#' genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
#' A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
#' B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
#' 
#'  
#' pt2 <- permTest(A=A, B=B, ntimes=10, alternative="auto", verbose=TRUE, genome=genome, evaluate.function=meanDistance, randomize.function=randomizeRegions, non.overlapping=FALSE)
#' summary(pt2)
#' plot(pt2)
#' plot(pt2, plotType="Tailed")
#' 
#' @references 
#' Davison, A. C. and Hinkley, D. V. (1997) Bootstrap methods and their application, Cambridge University Press, United Kingdom, 156-160
#' 
#' @export permTest


#min.parallel is specifies the minimum amount of work (number of regions in A per ntimes) to activate parallelization.


permTest <- function(A, ntimes=100, randomize.function, evaluate.function, alternative="auto", min.parallel=1000, force.parallel=NULL, randomize.function.name=NULL, evaluate.function.name=NULL, verbose=FALSE, ...) {
  
  #check arguments
  alternative<-match.arg(alternative,c("less","greater", "auto"))
  if(!hasArg(A)) stop("A is missing")
  if(!is.numeric(ntimes)) stop("ntime must be numeric")
  if(!hasArg(randomize.function)) stop("randomize.function is missing")
  if(!is.function(randomize.function)) stop("randomize.function must be a function")
  if(!hasArg(evaluate.function)) stop("evaluate.function is missing")
  if(!is.function(evaluate.function)) stop("evaluate.function must be a function")
  if(!is.numeric(min.parallel)) stop("min.parallel must be numeric")
  if(ntimes<100) print(paste0("Note: The minimum p-value with only ",ntimes," permutations is ",1/(ntimes+1),". You should consider increasing the number of permutations."))
  
  A <- toGRanges(A) #does nothing if already a GRanges object

 
  if(!is.null(force.parallel)) {
    doParallel <- force.parallel
  } else {
    doParallel <- (length(A)*ntimes > min.parallel)
  }
  
   
  #get the function names
  if(is.null(randomize.function.name)) {
    randomize.function.name <- match.call()["randomize.function"]
  }   
  if(is.null(evaluate.function.name)) {
    evaluate.function.name <- match.call()["evaluate.function"]
  }  
  
  
  
  #evaluate the A region set
  original.evaluate <- evaluate.function(A,...)
  
  if(!is.numeric(original.evaluate)) {
    stop(paste0("The evaluation function must return a numeric value but it returned an object of class ", class(original.evaluate)))
  }
  
  if(verbose) {
    #WARNING: to give some visual information about the computation done, we use a Progress Bar. However, since the GUI will not be updated
    #whilst in multi core processing (mclapply), we partition the work in chunks and repeatedly call the parallel computation, once per chunk.
    #Between chunks, the progess bar will be updated. The problem is that the total computation time might increase due to waiting and joining
    #multiple times.
    #Create the progress bar
  
    pb <- txtProgressBar(min = 0, max = ntimes, style = 3)
    setTxtProgressBar(pb, 0)
  }
 
  #define the function to create and evaluate the random sets
  randomize_and_evaluate <- function(foo, ...) {
    #randomize
    randomA <- randomize.function(A,...)
    #evaluate the random region set    
    if(verbose) {
      setTxtProgressBar(pb, foo)
    }
    
    return(evaluate.function(randomA,...))
  }
  
  #create the random sets and evaluate them  
  if(doParallel) {
    if(verbose) { #if verbose, we will do the computations in chunks and update the progress bar in between
      random.evaluate <- numeric()
      chunk.size <- max(round(ntimes/100+1), 10)
      e <- 0
      done <- FALSE
      while(!done) {
        s <- e + 1
        e <- s + chunk.size
        if(e >= ntimes) {
          e <- ntimes
          done <- TRUE
        }
        random.evaluate <- c(random.evaluate, unlist(mclapply(c(s:e), randomize_and_evaluate, ...)))
        setTxtProgressBar(pb, e)
      }    
    } else { #if not verbose, just do it
      random.evaluate <- unlist(mclapply(c(1:ntimes), randomize_and_evaluate, ...))
    }
  } else {
    random.evaluate <- unlist(lapply(c(1:ntimes), randomize_and_evaluate, ...))
  }
  
  
  nas<-length(which(is.na(random.evaluate)))
  if(nas>0) warning(paste0(nas," iterations returned NA's. Only ",ntimes-nas," iterations have been used to compute the p-value."))
  
  
  if(alternative == "auto") {
    if(original.evaluate < mean(random.evaluate,na.rm=TRUE)) {
      alternative <- "less"
    } else {
      alternative <- "greater"
    }
  }
  
  #Perform the statistical test and return the values generated 
  if (alternative == "less")    pval <- (sum(original.evaluate > random.evaluate, na.rm=TRUE) + 1) / (ntimes - nas + 1)
  if (alternative == "greater")    pval <- (sum(original.evaluate < random.evaluate, na.rm=TRUE) + 1) / (ntimes - nas + 1)
  
  if(original.evaluate == 0 & all(random.evaluate == 0)){
    warning(paste0("All permuted values and the original evaluation value are equal to 0. Z-score cannot be computed."))
    pval <- 1
    zscore <- NA
  } else{
    zscore <- round((original.evaluate - mean(random.evaluate, na.rm=TRUE)) / sd(random.evaluate, na.rm=TRUE), 4)
  }
  
    
  #Create the return object
  res<-list(pval=pval, ntimes=ntimes, alternative=alternative, observed=original.evaluate, permuted=random.evaluate, zscore=zscore,
            evaluate.function=evaluate.function, evaluate.function.name=evaluate.function.name,
            randomize.function=randomize.function, randomize.function.name=randomize.function.name)
  
  
  if(!is.finite(zscore)){
    warning(paste0("All permuted values are equal to ",random.evaluate[1],". Z-score is infinite."))
  }      
  

  
  if(alternative=="greater" & original.evaluate<mean(random.evaluate,na.rm=TRUE)) warning("Alternative is greater and the observed statistic is less than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
  if(alternative=="less" & original.evaluate>mean(random.evaluate,na.rm=TRUE)) warning("Alternative is less and the observed statistic is greater than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
  
  
  class(res) <- "permTestResults"
  return(res)

}
