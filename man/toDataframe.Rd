% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toDataframe.R
\name{toDataframe}
\alias{toDataframe}
\title{toDataframe}
\usage{
toDataframe(A, stranded=FALSE)
}
\arguments{
\item{A}{a \code{\link{GRanges}} object.}

\item{stranded}{(only used when A is a \code{\link{GRanges}} object) a logical indicating whether a column with the strand information have to be added to the result (Defaults to FALSE)}
}
\value{
A \code{data.frame} with the regions in A. If A was a \code{\link{GRanges}} object, the output will include any metadata present in A.
}
\description{
Transforms a \code{\link{GRanges}} object or a \code{\link{data.frame}}containing a region set into a \code{\link{data.frame}}.
}
\details{
If the oject is of class \code{\link{data.frame}}, it will be returned untouched.
}
\examples{
A <- data.frame(chr=1, start=c(1, 15, 24), end=c(10, 20, 30), x=c(1,2,3), y=c("a", "b", "c"))

A2 <- toGRanges(A)

toDataframe(A2)

}
\seealso{
\code{\link{toGRanges}}
}
