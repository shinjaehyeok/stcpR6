test_that("Less and two.sided - 1. Normal", {
  mu <- 1
  dl <- 0.1
  du <- 3
  # Initialization
  stcp_greater <- Stcp$new(
    method = "ST",
    family = "Normal",
    alternative = "greater",
    threshold = log(1 / 0.025),
    m_pre = mu,
    delta_lower = dl,
    delta_upper = du
  )

  stcp_less <- Stcp$new(
    method = "ST",
    family = "Normal",
    alternative = "less",
    threshold = log(1 / 0.025),
    m_pre = mu,
    delta_lower = dl,
    delta_upper = du
  )
  
  stcp_two_sided <- Stcp$new(
    method = "ST",
    family = "Normal",
    alternative = "two.sided",
    threshold = log(1 / 0.05),
    m_pre = mu,
    delta_lower = dl,
    delta_upper = du
  )
  
  # obs
  xs <- c(-1, 0, 2, 1)
  
  updates_greater <- stcp_greater$updateAndReturnHistories(xs)
  updates_less <- stcp_less$updateAndReturnHistories(xs)
  
  update_two_sided <- stcp_two_sided$updateAndReturnHistories(xs)

  # Two_sided is the average of greater and less updates
  expected_two_sided <- log(exp(updates_greater) + exp(updates_less)) - log(2)
  expect_equal(update_two_sided, expected_two_sided)  
  
  # less updates are equal to greater updates with flipped observaions.
  # flipped obs
  ys <- -xs
  stcp_greater_flipped <- Stcp$new(
    method = "ST",
    family = "Normal",
    alternative = "greater",
    threshold = log(1 / 0.025),
    m_pre = -mu,
    delta_lower = dl,
    delta_upper = du
  )
  
  expected_less <- stcp_greater_flipped$updateAndReturnHistories(ys)
  expect_equal(updates_less, expected_less)  
})

test_that("Less and two.sided - 2. Ber", {
  p <- 0.3
  dl <- 0.1
  du <- 0.2
  # Initialization
  stcp_greater <- Stcp$new(
    method = "ST",
    family = "Ber",
    alternative = "greater",
    threshold = log(1 / 0.025),
    m_pre = p,
    delta_lower = dl,
    delta_upper = du
  )
  
  stcp_less <- Stcp$new(
    method = "ST",
    family = "Ber",
    alternative = "less",
    threshold = log(1 / 0.025),
    m_pre = p,
    delta_lower = dl,
    delta_upper = du
  )
  
  stcp_two_sided <- Stcp$new(
    method = "ST",
    family = "Ber",
    alternative = "two.sided",
    threshold = log(1 / 0.05),
    m_pre = p,
    delta_lower = dl,
    delta_upper = du
  )
  
  # obs
  xs <- c(0, 1, 0, 0)
  
  updates_greater <- stcp_greater$updateAndReturnHistories(xs)
  updates_less <- stcp_less$updateAndReturnHistories(xs)
  
  update_two_sided <- stcp_two_sided$updateAndReturnHistories(xs)
  
  # Two_sided is the average of greater and less updates
  expected_two_sided <- log(exp(updates_greater) + exp(updates_less)) - log(2)
  expect_equal(update_two_sided, expected_two_sided)  
  
  # less updates are equal to greater updates with flipped observations.
  # flipped obs
  ys <- 1-xs
  stcp_greater_flipped <- Stcp$new(
    method = "ST",
    family = "Ber",
    alternative = "greater",
    threshold = log(1 / 0.025),
    m_pre = 1-p,
    delta_lower = dl,
    delta_upper = du
  )
  expected_less <- stcp_greater_flipped$updateAndReturnHistories(ys)
  expect_equal(updates_less, expected_less)  
})

test_that("Less and two.sided - 3. Bounded", {
  p <- 0.3
  dl <- 0.1
  du <- 0.2
  # Initialization
  stcp_greater <- Stcp$new(
    method = "ST",
    family = "Bounded",
    alternative = "greater",
    threshold = log(1 / 0.025),
    m_pre = p,
    delta_lower = dl,
    delta_upper = du
  )
  
  stcp_less <- Stcp$new(
    method = "ST",
    family = "Bounded",
    alternative = "less",
    threshold = log(1 / 0.025),
    m_pre = p,
    delta_lower = dl,
    delta_upper = du
  )
  
  stcp_two_sided <- Stcp$new(
    method = "ST",
    family = "Bounded",
    alternative = "two.sided",
    threshold = log(1 / 0.05),
    m_pre = p,
    delta_lower = dl,
    delta_upper = du
  )
  
  # obs
  xs <- c(0, 1, 0, 0)
  
  updates_greater <- stcp_greater$updateAndReturnHistories(xs)
  updates_less <- stcp_less$updateAndReturnHistories(xs)
  
  update_two_sided <- stcp_two_sided$updateAndReturnHistories(xs)
  
  # Two_sided is the average of greater and less updates
  expected_two_sided <- log(exp(updates_greater) + exp(updates_less)) - log(2)
  expect_equal(update_two_sided, expected_two_sided)  
  
  # less updates are equal to greater updates with flipped observations.
  # flipped obs
  ys <- 1-xs
  stcp_greater_flipped <- Stcp$new(
    method = "ST",
    family = "Bounded",
    alternative = "greater",
    threshold = log(1 / 0.025),
    m_pre = 1-p,
    delta_lower = dl,
    delta_upper = du
  )
  expected_less <- stcp_greater_flipped$updateAndReturnHistories(ys)
  expect_equal(updates_less, expected_less)  
})
