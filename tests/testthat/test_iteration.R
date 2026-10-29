test_that("permutation evaluation supports integer and double statistics", {
  A <- toGRanges("chr1:100-110")
  state <- new.env(parent=emptyenv())
  state$iteration <- 0L
  randomize <- function(A, offset) {
    state$iteration <- state$iteration + 1L
    GenomicRanges::shift(A, state$iteration * offset)
  }
  evaluate <- list(
    integer=function(A, offset) as.integer(sum(start(A)) + offset),
    double=function(A, offset) sum(start(A)) / 2 + offset
  )
  result <- permTest(A, ntimes=100, randomize.function=randomize,
                     evaluate.function=evaluate, force.parallel=FALSE,
                     offset=2)

  expect_named(result, c("integer", "double"))
  expect_equal(result$integer$observed, 102)
  expect_equal(result$integer$permuted, 102 + 2 * seq_len(100))
  expect_equal(result$double$observed, 52)
  expect_equal(result$double$permuted, 52 + seq_len(100))
})

test_that("a local z-score window smaller than its step evaluates only zero", {
  A <- toGRanges("chr1:100-110")
  pt <- structure(list(permuted=seq_len(5), zscore=0,
                       evaluate.function=length), class="permTestResults")
  result <- localZScore(A, pt, window=1, step=5)

  expect_equal(result$shifts, 0)
  expect_length(result$shifted.z.scores, 1)
})

test_that("resampling supports empty and single-region universes", {
  A <- toGRanges("chr1:100-110")
  empty <- A[integer()]

  expect_equal(resampleRegions(A, A), A)
  expect_equal(resampleRegions(A, A, per.chromosome=TRUE), A)
  expect_length(resampleRegions(empty, empty), 0)
  expect_length(resampleRegions(empty, empty, per.chromosome=TRUE), 0)
})
