test_that("Confidence sequence for Normal - 1. Simple", {
  # If lambda = 1 and sig = 1 then the corresponding CS is given by
  # x_bar - 1/2 - log(1/alpha) / n
  
  alpha <- 0.05
  n <- 100
  x_bar <- 0
  expected_width <- 1/2 + log(1/alpha) / n
  expected_interval <- c(x_bar - expected_width, Inf)

  # Note lambda = 1 corresponds to n = 2 * log(1/alpha)
  # from the relationship log(1/alpha) / n = psi_star(delta) = 1/2 for
  # the simple normal case.
  
  # Greater
  cs_greater <- NormalCS$new(
    alternative = "greater",
    alpha = alpha,
    n_upper = 2 * log(1/alpha),
    n_lower = 2 * log(1/alpha),
    weights = NULL,
    lambdas = NULL,
    k_max = 1000
  )
  
  expect_equal(cs_greater$getLambdas(), 1)
  expect_equal(cs_greater$computeWidth(n), expected_width)
  expect_equal(cs_greater$computeWidth(n, sig = 2), 2 * expected_width)
  
  cs_interval <- cs_greater$computeInterval(n, sig = 2)
  expect_equal(cs_interval, c(x_bar - 2 * expected_width, Inf))
  
  # Less
  cs_less <- NormalCS$new(
    alternative = "less",
    alpha = alpha,
    n_upper = 2 * log(1/alpha),
    n_lower = 2 * log(1/alpha),
    weights = NULL,
    lambdas = NULL,
    k_max = 1000
  )
  
  expect_equal(cs_less$getLambdas(), -1)
  expect_equal(cs_less$computeWidth(n), expected_width)
  expect_equal(cs_less$computeWidth(n, sig = 2), 2 * expected_width)
  
  cs_interval <- cs_less$computeInterval(n, sig = 2)
  expect_equal(cs_interval, c(-Inf, x_bar + 2 * expected_width))
  
  # Two.sided
  cs_two <- NormalCS$new(
    alternative = "two.sided",
    alpha = alpha,
    n_upper = 2 * log(1/alpha),
    n_lower = 2 * log(1/alpha),
    weights = NULL,
    lambdas = NULL,
    k_max = 1000
  )
  
  upper_width <- 1/2 + log(2/alpha) / n # Nixture gets a slightly narrower width
  
  expect_equal(cs_two$getLambdas(), c(1, -1))
  expect_true(cs_two$computeWidth(n) <= upper_width)
  expect_true(cs_two$computeWidth(n, sig = 2) <= 2 * upper_width)

  cs_interval <- cs_two$computeInterval(n, sig = 2)
  expect_true(sum(abs(cs_interval - c(x_bar - 2 * upper_width, x_bar + 2 * upper_width))) < 0.01)
  
})