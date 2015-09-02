library(regioneR)
context("Built-in Randomization Functions")

#Define some GRanges to use in the tests
emptyGR <- toGRanges(data.frame(chr=character(), start=numeric(), end=numeric()))
smallA <- toGRanges(data.frame(chr=rep(c("chr1", "chr2"), 10), start=100*(1:20), end=120*(1:20)))
smallB <- toGRanges(data.frame(chr=rep(c("chr2", "chr1"), 10), start=100*(1:20), end=105*(1:20)))

bigRegionsA <- toGRanges(data.frame(chr=rep(c("chr1", "chr2"), 10), start=100*(1:20), end=1200000*(1:20)))

gam <- getGenomeAndMask("hg19")

#Randomize Regions
test_that("the class of randomized regions is correct (randomizeRegions)", {
  expect_is(randomizeRegions(smallA), "GRanges")
  expect_is(randomizeRegions(emptyGR), "GRanges")
})

test_that("the number of randomized regions is correct (randomizeRegions)", {
  expect_equal(length(randomizeRegions(smallA)), length(smallA))
  expect_equal(length(randomizeRegions(smallA, per.chromosome=TRUE)), length(smallA))
  expect_equal(length(randomizeRegions(smallA, non.overlapping=TRUE)), length(smallA))

})

test_that("the randomized regions do not overlap the mask (randomizeRegions)", {
  expect_equal(numOverlaps(randomizeRegions(bigRegionsA, genome=gam$genome, mask=gam$mask), gam$mask), 0)
  expect_equal(numOverlaps(randomizeRegions(bigRegionsA, genome=gam$genome, mask=gam$mask, per.chromosome=TRUE), gam$mask), 0)
  expect_equal(numOverlaps(randomizeRegions(bigRegionsA, genome=gam$genome, mask=gam$mask, non.overlapping=TRUE), gam$mask), 0)    
})
          
test_that("the randomized regions do not overlap betwen them when allow.overlaps=FALSE (randomizeRegions)", {
  overlapsItself <- function(A) {
    ff <- findOverlaps(A, A)
    overlapping.regs <- unique(queryHits(ff)[((queryHits(ff) - subjectHits(ff)) != 0)]) #which region overlaps any region that is not itself?
    return(length(overlapping.regs)>0)
  }
  expect_false(overlapsItself(randomizeRegions(bigRegionsA, allow.overlaps=FALSE)))
  expect_false(overlapsItself(randomizeRegions(bigRegionsA[1:10], per.chromosome=TRUE, allow.overlaps=FALSE)))
})
          
     