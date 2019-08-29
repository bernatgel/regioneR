#' Permutation Test
#' 
#' Performs a permutation test to see if there is an association between a region set and some other feature using
#' an evaluation function.
#'  
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
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom methods hasArg is
#' @importFrom S4Vectors subjectHits queryHits
#' 
#' 
#' @export permTest
#' 
#' @importFrom stats sd


#min.parallel is specifies the minimum amount of work (number of regions in A per ntimes) to activate parallelization.


permTest <- function(A, ntimes=100, randomize.function, evaluate.function, alternative="auto", min.parallel=1000, force.parallel=NULL, randomize.function.name=NULL, evaluate.function.name=NULL, verbose=FALSE, ...) {
  
  #check arguments
  alternative<-match.arg(alternative,c("less","greater", "auto"))
  if(!hasArg(A)) stop("A is missing")
  if(!is.numeric(ntimes)) stop("ntime must be numeric")
  if(!hasArg(randomize.function)) stop("randomize.function is missing")
  if(!is.function(randomize.function)) stop("randomize.function must be a function")
  if(!hasArg(evaluate.function)) stop("evaluate.function is missing")
  if(!(is.function(evaluate.function) | is.list(evaluate.function))) stop("evaluate.function must be a function")
  if(!is.numeric(min.parallel)) stop("min.parallel must be numeric")
  if(ntimes<100) print(paste0("Note: The minimum p-value with only ",ntimes," permutations is ",1/(ntimes+1),". You should consider increasing the number of permutations."))
  
  A <- toGRanges(A) #does nothing if already a GRanges object

 
  if(!is.null(force.parallel)) {
    doParallel <- force.parallel
  } else {
    doParallel <- (length(A)*ntimes > min.parallel)
    if(verbose) message("Auto-setting parallel computing to: ", doParallel)
  }
  
    
  #Evaluation Function: get the function name and convert to list if its not yet
  if(!is.list(evaluate.function)) { #if it's a single function
    if(is.null(evaluate.function.name)) {
      evaluate.function.name <- as.character(match.call()["evaluate.function"])
    }
    ef <- list()
    ef[[evaluate.function.name]] <-  evaluate.function
    evaluate.function <- ef
  } else { #if it's a list of functions
    if(!is.null(evaluate.function.name)) { #if names were explicitely provided
      names(evaluate.function) <- evaluate.function.name
    } else { #try to leave the current names or create new ones if no names present
      if(is.null(names(evaluate.function))) {
        names(evaluate.function) <- paste0("Function", seq_len(length(evaluate.function)))
      }
    }
  }
  
  #Randomization Function: Get a name
  if(is.null(randomize.function.name)) {
    randomize.function.name <- match.call()["randomize.function"]
  }
      
  #Start the permutation test
  #compute the evaluation function(s) using the original region set A
  original.evaluate <- sapply(seq_len(length(evaluate.function)), function(i,...) {return(evaluate.function[[i]](A,...))}, ...)
 
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
    
    #compute the evaluation function(s) using the RANDOMIZED region set randomA
    rand.evaluate <- sapply(seq_len(length(evaluate.function)), function(i, ...) {return(evaluate.function[[i]](randomA,...))}, ...)
    
    return(rand.evaluate)
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
        random.evaluate <- rbind(random.evaluate, do.call(rbind, mclapply(c(s:e), randomize_and_evaluate, ...)))
        setTxtProgressBar(pb, e)
      }    
    } else { #if not verbose, just do it
      random.evaluate <- do.call(rbind, mclapply(seq_len(ntimes), randomize_and_evaluate, ...))
    }
  } else {
    random.evaluate <- do.call(rbind, lapply(seq_len(ntimes), randomize_and_evaluate, ...))
  }
  
 
  #The simulation process has finished. Now build a permTestResults object for each evaluate.function
  results <- list()
  for(i in seq_len(length(evaluate.function))) {
    #Get the data for the i-th function
    func.name <- names(evaluate.function)[i]
    orig.ev <- original.evaluate[i]
    rand.ev <- random.evaluate[,i]
    
    
    #warn if any NA or NaN
    num.nas <- length(which(is.na(rand.ev)))
    num.valid.values <- ntimes-num.nas
        
    if(num.valid.values < ntimes) {
      if(num.valid.values>0) {
        warning(paste0(num.nas, " iterations returned NA or NaN. Only ",  ," iterations have been used to compute the p-value."))
      } else {
        warning(paste0("All ", num.nas, " iterations returned NA or NaN. No valid values returned. It is not possible to compute the p-value nor z-score."))
      }
    }
    

    if(num.valid.values > 0) {
      #decide the alternative if alternative == "auto"
      if(alternative == "auto") {
        alt <- ifelse(orig.ev < mean(rand.ev, na.rm=TRUE), "less", "greater")
      } else {
        alt <- alternative
      }
            
      #Compute the p-value
      if (alt == "less") {
        pval <- (sum(orig.ev >= rand.ev, na.rm=TRUE) + 1) / (num.valid.values + 1)
      } else { #alt == "greater"
        pval <- (sum(orig.ev <= rand.ev, na.rm=TRUE) + 1) / (num.valid.values + 1)
      }
      #if the original alternative was not the best one, suggest the user to change it
      if(alternative=="greater" & orig.ev<mean(rand.ev,na.rm=TRUE)) message("Alternative is greater and the observed statistic is less than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
      if(alternative=="less" & orig.ev>mean(rand.ev,na.rm=TRUE)) message("Alternative is less and the observed statistic is greater than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
      #Compute the z-score
      if(orig.ev == 0 & all(rand.ev == 0)){ #If everything is 0, warning and "empty" results
        warning(paste0("All permuted values and the original evaluation value are equal to 0. Z-score cannot be computed."))
        pval <- 1
        zscore <- NA
      } else{
        zscore <- round((orig.ev - mean(rand.ev, na.rm=TRUE)) / stats::sd(rand.ev, na.rm=TRUE), 4)
      }
    } else {
      pval <- NA
      zscore <- NA
      alt <- alternative
    }
    
    
    if(!is.na(pval)) {
      if(!is.finite(zscore)){ #if all evaluations are equal, the sd is 0 and the z-score is infinite
        warning(paste0("All permuted values are equal to ", rand.ev[1], ". Z-score is infinite."))
      }  
    }
      
    #Create the permTestResults object
    res<-list(pval=pval, ntimes=ntimes, alternative=alt, observed=orig.ev, permuted=rand.ev, zscore=zscore,
                evaluate.function=evaluate.function[[i]], evaluate.function.name=func.name,
                randomize.function=randomize.function, randomize.function.name=randomize.function.name)
    class(res) <- "permTestResults"  
    results[[func.name]] <- res
  }

  class(results) <- "permTestResultsList"

  return(results)

}
