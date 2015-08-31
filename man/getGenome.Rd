\name{getGenome}
\alias{getGenome}
\title{getGenome}
\usage{
getGenome(...)
}
\arguments{
  \item{...}{a genome parameter is required. The genome
  object or genome identifier.}
}
\value{
A GRanges object with the "genome" data c(Chromosome, Start
(by default, 1), Chromosome Length) given a
\code{\link{BSgenome}}, a genome name, a
\code{\link{data.frame}} or a GRanges.

A \code{\link{GRanges}} representing the genome with one
region per chromosome.
}
\description{
Function to obtain a genome
}
\details{
If genome is a \code{\link{BSgenome}} (from the package
\code{BioStrings}), it will transform it into a
\code{\link{GRanges}} with chromosomes and chromosome
lengths.

If genome is a \code{\link{data.frame}} with 3 columns, it
will transform it into a GRanges.

If genome is a \code{\link{data.frame}} with 2 columns, it
will assume the first is the chromosome, the second is the
length of the chromosomes and will add 1 as start.

If genome is a \code{character} string uniquely identifying
a \code{\link{BSgenome}} installed in the system (e.g.
"hg19", "mm10",... but not "hg"), it will create a genome
based on the \code{\link{BSgenome}} object identified by
the character string.

If genome is a \code{ \link{GRanges}} object, it will
return it as is.

If genome is non of the above, it will give a warning and
try to transform it into a GRanges using \link{toGRanges}.
This can be helpful if \code{genome} is a connection to a
file.
}
\note{
This function is memoised (cached) using the
\code{\link{memoise}} package. To empty the cache, use
\code{\link{forget}(getGenome)}

Please note that passing this function the path to a file
will not work, since it will assume the character is the
identifier of a genome. To read the genome from a file,
please use \code{getGenome(toGRanges("path/to/file"))}
}
\examples{
getGenome("hg19")

getGenome(data.frame(c("chrA", "chrB"), c(15000000, 10000000)))
}
\seealso{
\code{\link{getMask}}, \code{\link{getGenomeAndMask}},
\code{\link{characterToBSGenome}},
\code{\link{maskFromBSGenome}},
\code{\link{emptyCacheRegioneR}}
}

