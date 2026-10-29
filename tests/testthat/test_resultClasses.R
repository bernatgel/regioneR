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
