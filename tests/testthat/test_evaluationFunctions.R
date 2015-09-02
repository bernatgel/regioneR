library(regioneR)
context("Built-in Evaluation Functions")

#Define some GRanges to use in the tests
emptyGR <- toGRanges(data.frame(chr=character(), start=numeric(), end=numeric()))
smallA <- toGRanges(data.frame(chr=rep(c("chr1", "chr2"), 10), start=100*(1:20), end=120*(1:20)))
smallB <- toGRanges(data.frame(chr=rep(c("chr2", "chr1"), 10), start=100*(1:20), end=105*(1:20)))

bigA <- toGRanges(system.file("extdata", "my.special.genes.bed", package="regioneR"))
universeA <- toGRanges(system.file("extdata", "all.genes.bed", package="regioneR"))
bigB <- toGRanges(system.file("extdata", "my.altered.regions.bed", package="regioneR"))

#Test numOverlaps
test_that("the numOverlaps function returns a correct result", {
  expect_equal(numOverlaps(smallA, smallA, count.once=FALSE), 38)
  expect_equal(numOverlaps(smallA, smallA, count.once=TRUE), length(smallA))
  expect_equal(numOverlaps(smallA, smallB), 18)  
  expect_equal(numOverlaps(smallA, smallB, count.once=TRUE), 15)  
  expect_equal(numOverlaps(emptyGR, emptyGR), 0)
  expect_equal(numOverlaps(smallA, emptyGR, count.once=FALSE), 0)  
  expect_equal(numOverlaps(smallA, emptyGR, count.once=TRUE), 0)  
  
  expect_equal(numOverlaps(bigA, universeA, count.once=TRUE), length(bigA))
  expect_equal(numOverlaps(universeA, bigA, count.once=FALSE), numOverlaps(universeA, bigA, count.once=F))
  expect_equal(numOverlaps(universeA, bigA, count.once=FALSE), numOverlaps(bigA, universeA))

  expect_equal(numOverlaps(smallA, emptyGR, foo=TRUE, bar=3, baz="C"), 0)
})

#meanDistance
test_that("the meanDistance function returns a correct result", {
  expect_equal(meanDistance(smallA, smallB), 10)
  expect_equal(meanDistance(smallA, smallA), 0)
  expect_equal(meanDistance(A=smallA, B=smallA), 0)
})
