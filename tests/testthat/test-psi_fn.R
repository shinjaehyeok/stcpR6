test_that("Sub-G psi functions", {
  psi_fn_list <- generate_sub_G_fn(1)
  
  # At zero
  expect_equal(psi_fn_list$psi(0), 0)
  expect_equal(psi_fn_list$psi_star(0), 0)
  expect_equal(psi_fn_list$psi_star_div(0), 0)
  expect_equal(psi_fn_list$psi_star_inv(0), 0)
  
  # Derivatives always positive except at zero
  y_vec <- seq(0.1, 10, by = 0.1)
  div_vec <- sapply(y_vec, psi_fn_list$psi_star_div)
  expect_equal(sum(div_vec > 0), length(y_vec))
  
  # Inverse function
  mu_vec <- seq(0, 10, by = 0.1)
  psi_star_vec <- sapply(mu_vec, psi_fn_list$psi_star)
  psi_inv_from_star <- sapply(psi_star_vec, psi_fn_list$psi_star_inv)
  expect_equal(sum(abs(mu_vec - psi_inv_from_star) > 1e-6), 0)
  
  lambda_vec <- seq(0, 10, by = 0.1)
  psi_star_inv_vec <- sapply(lambda_vec, psi_fn_list$psi_star_inv)  
  psi_star_from_inv <- sapply(psi_star_inv_vec, psi_fn_list$psi_star)
  expect_equal(sum(abs(lambda_vec - psi_star_from_inv) > 1e-6), 0)
})

test_that("Sub-B psi functions", {
  m_pre <- 0.6
  psi_fn_list <- generate_sub_B_fn(m_pre)
  
  # At zero
  expect_equal(psi_fn_list$psi(0), 0)
  expect_equal(psi_fn_list$psi_star(0), 0)
  expect_equal(psi_fn_list$psi_star_div(0), 0)
  expect_equal(psi_fn_list$psi_star_inv(0), 0)
  
  # Derivatives always positive except at zero
  y_vec <- seq(0.1, 1-m_pre, by = 0.01)
  div_vec <- sapply(y_vec, psi_fn_list$psi_star_div)
  expect_equal(sum(div_vec > 0), length(y_vec))
  
  # Inverse function
  mu_vec <- seq(0, 1-m_pre, by = 0.01)
  psi_star_vec <- sapply(mu_vec, psi_fn_list$psi_star)
  psi_inv_from_star <- sapply(psi_star_vec, psi_fn_list$psi_star_inv)
  expect_equal(sum(abs(mu_vec - psi_inv_from_star) > 1e-6), 0)

  lambda_vec <- seq(0, log(1/m_pre), by = 0.01)
  psi_star_inv_vec <- sapply(lambda_vec, psi_fn_list$psi_star_inv)
  psi_star_from_inv <- sapply(psi_star_inv_vec, psi_fn_list$psi_star)
  expect_equal(sum(abs(lambda_vec - psi_star_from_inv) > 1e-6), 0)
})

test_that("Sub-E psi functions", {
  psi_fn_list <- generate_sub_E_fn()
  
  # At zero
  expect_equal(psi_fn_list$psi(0), 0)
  expect_equal(psi_fn_list$psi_star(0), 0)
  expect_equal(psi_fn_list$psi_star_div(0), 0)
  expect_equal(psi_fn_list$psi_star_inv(0), 0)
  
  # Derivatives always positive except at zero
  y_vec <- seq(0.1, 10, by = 0.1)
  div_vec <- sapply(y_vec, psi_fn_list$psi_star_div)
  expect_equal(sum(div_vec > 0), length(y_vec))
  
  # Inverse function
  mu_vec <- seq(0, 10, by = 0.1)
  psi_star_vec <- sapply(mu_vec, psi_fn_list$psi_star)
  psi_inv_from_star <- sapply(psi_star_vec, psi_fn_list$psi_star_inv)
  expect_equal(sum(abs(mu_vec - psi_inv_from_star) > 1e-6), 0)
  
  lambda_vec <- seq(0, 10, by = 0.1)
  psi_star_inv_vec <- sapply(lambda_vec, psi_fn_list$psi_star_inv)  
  psi_star_from_inv <- sapply(psi_star_inv_vec, psi_fn_list$psi_star)
  expect_equal(sum(abs(lambda_vec - psi_star_from_inv) > 1e-6), 0)
})