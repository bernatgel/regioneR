\name{toGRanges}
\alias{toGRanges}
\title{toGRanges}
\usage{
toGRanges(A, ...)
}
\arguments{
  \item{A}{a \code{\link{data.frame}} containing a region
  set, a \code{\link{GRanges}} object, a BED file or any
  type of file supported by \code{rtracklayer}}

  \item{...}{further arguments to be passed to other
  methods.}
}
\value{
A \code{\link{GRanges}} object with the regions in A
}
\description{
Transforms a file or an object containing a region set into
a \code{\link{GRanges}} object.
}
\details{
If A is already a \code{\link{GRanges}} object, it will be
returned untouched.

If A is a file name or connection to a file in any of the
formats supported by \code{rtracklayer}'s import function
(BED, GFF...) it will be imported using \code{rtracklayer}.

If A is a data frame, the function will assume the first
three columns are chromosome, start and end and create a
\code{\link{GRanges}} object. Any additional column will be
considered metadata and stored as such in the
\code{\link{GRanges}} object.
}
\examples{
A <- data.frame(chr=1, start=c(1, 15, 24), end=c(10, 20, 30),  x=c(1,2,3), y=c("a", "b", "c"))

toGRanges(A)
}
\seealso{
\code{\link{toDataframe}}
}

