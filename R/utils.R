


getSeparator <- function(l, seps=c("\t", ",", ";", "|")) {
  #Check all seps and find the one where all rows are split into the same 
  #number of elements and said number is bigger than 1 and it's the largest
  num.fields <- unlist(lapply(seps, function(sep) {
    split.length <- unlist(lapply(strsplit(x = l, split = sep), length))
    if(split.length[1]>1 && all(split.length==split.length[1])) return(split.length[1])
    else return(-1)
  }))
  max.sep <- IRanges::which.max(num.fields)[1]  #if more than one sep is equally good, chose the first one
  if(num.fields[max.sep]>1) return(seps[max.sep])
  return(NULL)
}
 

hasHeader <- function(l, sep) {
  #Check if the first line does not contain numerics (since this is expected 
  #to be used by toGRanges, we expect there will be numeric columns)
  nums <- suppressWarnings(as.numeric(strsplit(x = l[1], split = sep)[[1]]))
  return(all(is.na(nums)))
}
  
  