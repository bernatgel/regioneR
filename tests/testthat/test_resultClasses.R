test_that("permutation result methods accept S3 subclasses", {
  pt <- structure(list(observed=4, permuted=1:5, ntimes=5,
                       alternative="greater", pval=0.5, zscore=0,
                       evaluate.function=length,
                       evaluate.function.name="length",
                       randomize.function.name="identity"),
                  class="permTestResults")
  derived <- pt
  class(derived) <- c("derivedPermTestResults", "permTestResults")
  A <- toGRanges("chr1:100-110")

  expect_equal(localZScore(A, derived, window=10, step=5),
               localZScore(A, pt, window=10, step=5))
  expected <- recomputePermTest(pt)
  class(expected) <- class(derived)
  expect_equal(recomputePermTest(derived), expected)
  expect_output(summary(derived), "Number of permutations: 5")

  results <- structure(list(derived),
                       class=c("derivedResultsList", "permTestResultsList"))
  expect_s3_class(localZScore(A, results, window=10, step=5),
                  "localZScoreResultsList")
  expect_output(summary(results), "Permutation tests: 1")
  empty <- structure(list(), class=class(results))
  expect_output(summary(empty),
                "permTestResultList object of length 0")

  grDevices::pdf(NULL)
  tryCatch(expect_no_error(plot(derived)), finally=grDevices::dev.off())
})

test_that("graphical overlap summary preserves its invisible option return", {
  A <- data.frame(chr=1, start=c(1, 15, 24, 40, 50), end=c(10, 20, 30, 45, 55))
  B <- data.frame(chr=1, start=c(2, 12, 28, 35), end=c(5, 25, 33, 43))
  old.scipen <- getOption("scipen")
  grDevices::pdf(NULL)
  result <- tryCatch(withVisible(overlapGraphicalSummary(A, B)),
                     finally=grDevices::dev.off())
  expect_false(result$visible)
  expect_equal(result$value, list(scipen=999))
  expect_identical(getOption("scipen"), old.scipen)
})

test_that("permutation result methods reject unrelated S3 classes", {
  unrelated <- structure(list(), class=c("unrelated", "list"))
  expect_error(localZScore(toGRanges("chr1:100-110"), unrelated),
               "pt must be of class permTestResults")
  expect_error(recomputePermTest(unrelated), "must be a permTestResults object")
  expect_error(regioneR:::summary.permTestResults(unrelated),
               "must be a permTestResults object")
  expect_error(regioneR:::summary.permTestResultsList(unrelated),
               "must be a permTestResultsList object")
})
