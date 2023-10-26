test_that("logSumExpTrick", {
  x <- 3
  expect_equal(logSumExpTrick(x), x)
  
  xs <- c(-1, 2, 0)
  expect_equal(logSumExpTrick(xs), log(sum(exp(xs))))
  
  expect_error(logSumExpTrick())
})

