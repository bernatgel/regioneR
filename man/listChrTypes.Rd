\name{listChrTypes}
\alias{listChrTypes}
\title{filterChromosomes
listChrTypes}
\usage{
listChrTypes()
}
\value{
the list of available chrs and organisms is printed
}
\description{
Prints a list of the available organisms and chromosomes
sets in the predefined chromosomes sets information.
}
\examples{
g <- getGenomeAndMask("hg19")$genome

listChrTypes()

g <- filterChromosomes(g, chr.type="autosomal", organism="hg19")
}
\seealso{
\code{\link{filterChromosomes}},
\code{\link{getChromosomesByOrganism}}
}

