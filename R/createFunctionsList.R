#' Create Functions List
#' 
#' @description 
#' Partially applies (the standard Curry function in functional programming) a list of arguments
#' to a function and returns a list of preapplied functions. The result of this function is a
#' list of functions suitable for the multiple evaluation functions in permTest.
#' 
#' 
#' @note 
#' It uses the code posted by "hadley" at http://stackoverflow.com/questions/6547219/how-to-bind-function-arguments
#' 
#' @usage 
#' createFunctionsList(FUN, param.name, values, func.names)
#' 
#' @param FUN   Function. the function to be partially applied
#' @param param.name  Character. The name of the parameter to pre-set.
#' @param values  A list or vector of values to preassign. A function will be created for each of the values in values. If present, the names of the list will be the names of the functions.
#' @param func.names  Character. The names of the functions created. Useful to identify the functions created. Defaults to the names of the values list or to Function1, Function2... if the values list has no names.
#' 
#' @return
#' It returns a list of functions with parameter param.value pre-set to values.
#' 
#' @seealso \code{\link{permTest}}, \code{\link{overlapPermTest}}
#' 
#' @examples
#' f <- function(a, b) {
#'  return(a+b)
#' }
#' 
#' funcs <- createFunctionsList(FUN=f, param.name="b", values=c(1,2,3), func.names=c("plusone", "plustwo", "plusthree"))
#' 
#' funcs$plusone(2)
#' funcs$plusone(10)
#' funcs$plusthree(2)
#' 
#' A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=0, mask=NA)
#' B <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=0, mask=NA)
#' 
#'  overlapsWith <- createFunctionsList(FUN=numOverlaps, param.name="B", values=list(a=A, b=B))
#'  overlapsWith$a(A=A) 
#'  overlapsWith$b(A=A) 
#' 
#' @export createFunctionsList


createFunctionsList <- function(FUN, param.name, values, func.names=NULL) {
  if(!hasArg(FUN)) stop("FUN is missing")
  if(!hasArg(param.name)) stop("param.name is missing")
  if(!is.character(param.name)) stop("param.name must be a character")
  if(length(param.name)>1) stop("param.name must be a single character string")
  if(!hasArg(values)) stop("values is missing")
  if(!(is.list(values) | is.vector(values))) stop("values must be a list or a vector")
  
  
  if(is.null(func.names)) {
    if(is.list(values) && !is.null(names(values))) {
      if(any(duplicated(names(values)))) stop("the names of the values list must be unique")
      func.names <- names(values)
    } else {
      func.names <- paste0("Function", c(1:length(values)))
    }
  }
  
  curried.funcs <- list()
  for(i in seq_len(length(values))) {
    curry.args <- list(FUN=FUN)
    curry.args[[param.name]] <- values[[i]]
    curried.funcs[[func.names[i]]] <- do.call("Curry", curry.args)    
  }
  return(curried.funcs)
}



#Curry function. From http://stackoverflow.com/questions/6547219/how-to-bind-function-arguments
Curry <- function(FUN, ...) {
  args <- match.call(expand.dots = FALSE)$...
  args$... <- as.name("...")
  
  env <- new.env(parent = parent.frame())
  
  if (is.name(FUN)) {
    fname <- FUN
  } else if (is.character(FUN)) {
    fname <- as.name(FUN)
  } else if (is.function(FUN)){
    fname <- as.name("FUN")
    env$FUN <- FUN
  } else {
    stop("FUN not function or name of function")
  }
  curry_call <- as.call(c(list(fname), args))
  
  f <- eval(call("function", as.pairlist(alist(... = )), curry_call))
  environment(f) <- env
  f
}