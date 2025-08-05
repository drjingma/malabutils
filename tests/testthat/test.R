# tests/testthat/test-rlog.R
test_that("rlog transforms non-zero values correctly", {
  x <- c(0, 1, 2, 4)
  expect_equal(rlog(x)[2], log(1))
  expect_equal(rlog(x)[3], log(2))
  expect_equal(rlog(x)[1], 0)
})
