test_that("GLRCU - Normal run as same as the hard-coded R function", {
  # Config
  mu_pre <- 1 # H0: mu = 1 H1: mu != 1
  max_sample <- 1000L
  ARL_target <- max_sample * 0.5
  
  # Hard-coded GLRCU Normal for "greater" alternative
  # We return the entire run history for a test
  glrcuNormalwithHistory <- function(x_vec, 
                                     mu_pre, 
                                     thres) {
    n_max <- length(x_vec)
    n_star <- Inf
    m <- -Inf
    mu_hat <- 0
    m_vec <- length(n_max)
    mu_hat_vec <- length(n_max)
    for (i in 1:n_max) {
      for (j in 1:i) {
        n_inner <- i - j + 1
        s <- sum(x_vec[j:i])
        mu_hat_inner <- max(s/n_inner, mu_pre)
        m_inner <-  n_inner * (mu_hat_inner - mu_pre)^2 / 2.0
        if (m < m_inner) {
          m <- m_inner
          mu_hat <- mu_hat_inner
        }
      }
      m_vec[i] <- m
      mu_hat_vec[i] <- mu_hat
      if (m > thres) {
        n_star <- i
      }
      m <- -Inf
    }
    return(list(n_star = n_star, m_vec = m_vec, mu_hat_vec = mu_hat_vec))
  }
  
  # stcpR6 implementation
  glrcuNormalStcp <- stcpR6::Stcp$new(
    method = "GLRCU",
    family = "Normal",
    alternative = "greater",
    threshold = log(ARL_target),
    m_pre = mu_pre
  )
  
  # Test sample
  x_vec <- rnorm(max_sample, mu_pre)
  
  # Runs
  glr_normal_R_out <- glrcuNormalwithHistory(x_vec, mu_pre, thres = log(ARL_target))
  
  glr_normal_Stcp_out <- glrcuNormalStcp$updateAndReturnHistories(x_vec)

  # Two runs are expected to be almost equal to each other
  # plot(glr_normal_R_out$m_vec, glr_normal_Stcp_out)
  expect_true(mean(abs(glr_normal_R_out$m_vec - glr_normal_Stcp_out)) < 1e-8)
})

test_that("GLRCU - Bernoulli run as same as the hard-coded R function", {
  # Config
  p_pre <- 0.3 # H0: p = 0.3 H1: p != 0.3
  max_sample <- 1000L
  ARL_target <- max_sample * 0.5

  # Hard-coded GLRCU Bernoulli for "two.sided" alternative
  # We return the entire run history for a test
  glrcuBerwithHistory <- function(x_vec, 
                                   p_pre, 
                                   thres) {
    n_max <- length(x_vec)
    n_star <- Inf
    m <- -Inf
    p_lower <- 1e-5
    p_upper <- 1 - p_lower
    p_hat <- 0.5
    m_vec <- length(n_max)
    p_hat_vec <- length(n_max)
    for (i in 1:n_max) {
      for (j in 1:i) {
        n_inner <- i - j + 1
        s <- sum(x_vec[j:i])
        f <- n_inner - s
        # For simplicity, we use a soft cap around boundary
        # which is slightly different from a formal Stcp implementation
        p_hat_inner <- min(max(s/n_inner, p_lower), p_upper)
        m_inner <-
          s * log(p_hat_inner/p_pre) + f * log((1-p_hat_inner)/(1-p_pre))
        if (m < m_inner) {
          m <- m_inner
          p_hat <- p_hat_inner
        }
      }
      m_vec[i] <- m
      p_hat_vec[i] <- p_hat
      if (m > thres) {
        n_star <- i
      }
      m <- -Inf
    }
    return(list(n_star = n_star, m_vec = m_vec, p_hat_vec = p_hat_vec))
  }
  
  # stcpR6 implementation
  glrcuBerStcp <- stcpR6::Stcp$new(
    method = "GLRCU",
    family = "Ber",
    alternative = "two.sided",
    threshold = log(ARL_target),
    m_pre = p_pre
  )
  
  # Test sample
  set.seed(1)
  x_vec <- rbinom(max_sample, 1, p_pre)
  
  # Runs
  glr_ber_R_out <- glrcuBerwithHistory(x_vec, p_pre, thres = log(ARL_target))
  
  glr_ber_Stcp_out <- glrcuBerStcp$updateAndReturnHistories(x_vec)

  # Two runs are expected to be almost equal to each other
  # Slight difference comes from p_hat implementation detail
  # plot(glr_ber_R_out$m_vec, glr_ber_Stcp_out)
  expect_true(mean(abs(glr_ber_R_out$m_vec - glr_ber_Stcp_out)) < 1e-4)
})