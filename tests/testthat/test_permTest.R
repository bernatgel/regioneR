library(regioneR)
context("Built-in Evaluation Functions")

#Define some GRanges to use in the tests
emptyGR <- toGRanges(data.frame(chr=character(), start=numeric(), end=numeric()))
smallA <- toGRanges(data.frame(chr=rep(c("chr1", "chr2"), 10), start=100*(1:20), end=120*(1:20)))
smallB <- toGRanges(data.frame(chr=rep(c("chr2", "chr1"), 10), start=100*(1:20), end=105*(1:20)))

bigA <- toGRanges(system.file("extdata", "my.special.genes.bed", package="regioneR"))
universeA <- toGRanges(system.file("extdata", "all.genes.bed", package="regioneR"))
bigB <- toGRanges(system.file("extdata", "my.altered.regions.bed", package="regioneR"))


#TODO
