test_that("Sequentail test for Normal - 1. Simple", {
  # Default: simple normal test N(0,1) vs N(1, 1) at level 0.05
  # Initialization
  stcp <- Stcp$new(
    method = "ST",
    family = "Normal",
    alternative = "greater",
    threshold = log(1 / 0.05),
    m_pre = 0,
    delta_lower = 1,
    delta_upper = 1
  )
  
  expect_equal(stcp$getWeights(), 1)
  expect_equal(stcp$getLambdas(), 1)
  
  expect_equal(stcp$getLogValue(), 0)
  expect_equal(stcp$getThreshold(), log(20))
  expect_equal(stcp$getTime(), 0)
  expect_equal(stcp$isStopped(), FALSE)
  expect_equal(stcp$getStoppedTime(), 0)
  
  
  obs <- c(1.0, 3.0, 2.0)
  expected_increment <- obs - 0.5 # x - 0.5
  expected_log_values <- cumsum(expected_increment)
  stcp$reset()
  updates <- stcp$updateAndReturnHistories(obs)
  expect_equal(updates, expected_log_values)
  expect_equal(stcp$getTime(), length(obs))
  expect_equal(stcp$isStopped(), TRUE)
  expect_equal(stcp$getStoppedTime(), 2)
  
  stcp$reset()
  stcp$updateLogValuesUntilStop(obs)
  expect_equal(updates, expected_log_values)
  expect_equal(stcp$getTime(), 2)
  expect_equal(stcp$isStopped(), TRUE)
  expect_equal(stcp$getStoppedTime(), 2)
})

test_that("Sequentail test for Normal - 2. Mixture", {
  # Initialization
  alpha <- 0.05
  weights <- c(0.3, 0.7)
  lambdas <- c(1, 2)
  stcp <- Stcp$new(
    method = "ST",
    family = "Normal",
    alternative = "greater",
    threshold = log(1 / alpha),
    m_pre = 0,
    weights = weights,
    lambdas = lambdas
  )
  
  expect_equal(stcp$getWeights(), weights)
  expect_equal(stcp$getLambdas(), lambdas)
  
  
  obs <- c(1.0, 3.0, 2.0)
  expected_increment1 <- lambdas[1] * obs - lambdas[1] ^ 2 / 2
  expected_increment2 <- lambdas[2] * obs - lambdas[2] ^ 2 / 2
  expected_log_values1 <- cumsum(expected_increment1)
  expected_log_values2 <- cumsum(expected_increment2)
  expected_log_value <- log(weights[1] * exp(expected_log_values1) +
                              weights[2] * exp(expected_log_values2))
  
  
  updates <- stcp$updateAndReturnHistories(obs)
  expect_equal(updates, expected_log_value)
  expect_equal(stcp$getTime(), length(obs))
  expect_equal(stcp$isStopped(), TRUE)
  expect_equal(stcp$getStoppedTime(), 2)
  
  stcp$reset()
  stcp$updateLogValuesUntilStop(obs)
  expect_equal(updates, expected_log_value)
  expect_equal(stcp$getTime(), 2)
  expect_equal(stcp$isStopped(), TRUE)
  expect_equal(stcp$getStoppedTime(), 2)
  
  
  x_bars <- c(2, 2)
  ns <- c(2, 1)
  expected_log_value_by_avgs <- expected_log_value[2:3]
  stcp$reset()
  updates <- stcp$updateAndReturnHistoriesByAvgs(x_bars, ns)
  expect_equal(updates, expected_log_value_by_avgs)
  expect_equal(stcp$getTime(), length(obs))
  expect_equal(stcp$isStopped(), TRUE)
  expect_equal(stcp$getStoppedTime(), 2)
  
})
