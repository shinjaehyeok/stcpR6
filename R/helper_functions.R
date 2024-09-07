#' log-sum-exp trick
#'
#' Apply log-sum-exp trick to a numeric vector.
#'
#' @param xs A numeric vector.
#'
#' @return log of sum of exp of \code{xs}, which is equal to \code{log(sum(exp(xs)))}.
#' @export
#'
logSumExpTrick <- function(xs) {
  return(logSumExp(xs))
}

#' Check whether the input delta parameters are acceptable
#'
#' For each method and family, check whether delta parameters are within
#' expected range with respect to the pre-change parameter.
#'
#' @param method Method of the sequential procedure.
#' * ST: Sequential test based on a mixture of E-values.
#' * SR: Sequential change detection based on e-SR procedure.
#' * CU: Sequential change detection based on e-CUSUM procedure.
#' * GLRCU: Sequential change detection based on GLR-CUSUM procedure.
#'
#' @param family Distribution of underlying univariate observations.
#' * Normal: (sub-)Gaussian with sigma = 1.
#' * Ber: Bernoulli distribution on \{0,1\}.
#' * Bounded: General bounded distribution on \[0,1\]
#'
#' @param alternative Alternative / post-change mean space
#' * two.sided: Two-sided test / change detection
#' * greater: Alternative /post-change mean is greater than null / pre-change one
#' * less:  Alternative /post-change mean is less than null / pre-change one
#'
#' @param m_pre The boundary of mean parameter in null / pre-change space
#'
#' @param delta_lower Minimum gap between null / pre-change space and
#' alternative / post-change one. It must be strictly positive for ST, SR and CU.
#' Currently, GLRCU does not support the minimum gap, and this param will be ignored.
#'
#' @param delta_upper Maximum gap between null / pre-change space and
#' alternative / post-change one. It must be strictly positive for ST, SR and CU.
#' Currently, GLRCU does not support the maximum gap, and this param will be ignored.
#'
#' @return A list of
#' 1. Boolean indicating whether it is acceptable or not.
#' 2. Character describing why it is not acceptable.
#' 3. Updated delta_upper for the case where the original input was NULL
#'
checkDeltaRange <- function(method,
                            family,
                            alternative,
                            m_pre,
                            delta_lower,
                            delta_upper) {
  is_acceptable <- FALSE
  error_message <- ""
  
  if (delta_lower <= 0 && method != "GLRCU") {
    stop("delta_lower must be positive.")
  }
  
  # Check delta_lower is within boundary for Ber and Bounded cases
  # For Ber, post-change parameter must be strictly within (0,1)
  # For Bounded case, post-change parameter can include 0 or 1.
  if (family == "Ber") {
    if (alternative != "less") {
      if (m_pre + delta_lower >= 1) {
        error_message <- paste0(
          "The minimum of alternative / post-change parameter ",
          "(m_pre + delta_lower) is greater than or equal to 1."
        )
      }
    }
    if (alternative != "greater") {
      if (m_pre - delta_lower <= 0) {
        error_message <- paste0(
          "The maximum of alternative / post-change parameter ",
          "(m_pre - delta_lower) is less than or equal to 0."
        )
      }
    }
  } else if (family == "Bounded") {
    if (alternative == "two.sided") {
      # Note if alternative == "greater",
      # bounded method technically does not require the upper bound 1.
      if (m_pre + delta_lower > 1) {
        error_message <- paste0(
          "The minimum of alternative / post-change parameter ",
          "(m_pre + delta_lower) is greater than 1."
        )
      }
    }
    if (alternative != "greater") {
      if (m_pre - delta_lower < 0) {
        error_message <- paste0(
          "The maximum of alternative / post-change parameter ",
          "(m_pre - delta_lower) is less than 0."
        )
      }
    }
  }
  
  # If delta_upper is NULL, we pick a reasonably large one
  if (is.null(delta_upper)) {
    if (family == "Normal") {
      delta_upper_ <- 5
    }
    if (family == "Ber" || family == "Bounded") {
      if (alternative == "greater") {
        delta_upper_ <- 1 - m_pre - 0.001
      } else if (alternative == "less") {
        delta_upper_ <- m_pre - 0.001
      } else {
        delta_upper_ <- min(1 - m_pre, m_pre) - 0.001
      }
    }
    delta_upper <- max(delta_lower, delta_upper_)
  }
  
  if (delta_upper <  delta_lower) {
    error_message <-
      "If not NULL, delta_upper must be greater than or equal to delta_lower."
  }
  
  is_acceptable <- ifelse(error_message == "", TRUE, FALSE)
  
  return(
    list(
      is_acceptable = is_acceptable,
      error_message = error_message,
      delta_upper = delta_upper
    )
  )
}

#' converted input deltas to parameters for exponential baselines
#'
#' For each exponential baseline family, convert delta range into
#' corresponding lambdas and weights.
#'
#' @param family Distribution of underlying univariate observations.
#' * Normal: (sub-)Gaussian with sigma = 1.
#' * Ber: Bernoulli distribution on \{0,1\}.
#' * Bounded: General bounded distribution on \[0,1\]
#'
#' @param alternative Alternative / post-change mean space
#' * two.sided: Two-sided test / change detection
#' * greater: Alternative /post-change mean is greater than null / pre-change one
#' * less:  Alternative /post-change mean is less than null / pre-change one
#'
#' @param threshold Stopping threshold. We recommend to use log(1/alpha)
#' for "ST" and "SR" methods where alpha is a testing level or 1/ARL.
#' for "CU" and "GRLCU", we recommend to tune the threshold by using
#' domain-specific sampler to hit the target ARL.
#'
#' @param m_pre The boundary of mean parameter in null / pre-change space
#'
#' @param delta_lower Minimum gap between null / pre-change space and
#' alternative / post-change one. It must be strictly positive for ST, SR and CU.
#' Currently, GLRCU does not support the minimum gap, and this param will be ignored.
#'
#' @param delta_upper Maximum gap between null / pre-change space and
#' alternative / post-change one. It must be strictly positive for ST, SR and CU.
#' Currently, GLRCU does not support the maximum gap, and this param will be ignored.
#'
#' @param k_max Positive integer to determine the maximum number of baselines.
#' For GLRCU method, it is used as the lookup window size for GLRCU statistics.
#'
#' @return A list of weights and lambdas
#' @export
#'
convertDeltaToExpParams <- function(family,
                                    alternative,
                                    threshold,
                                    m_pre,
                                    delta_lower,
                                    delta_upper,
                                    k_max) {
  alpha <- exp(-threshold)
  if (family == "Bounded") {
    # Bounded family uses sub-E class internally
    # So we convert it into the sub-E space.
    # where delta_E = m * delta / (sigma^2 + delta^2)
    delta_lower_internal_greater <-
      m_pre * delta_lower / (0.25 + delta_upper ^ 2)
    delta_upper_internal_greater <-
      m_pre * delta_upper / delta_lower ^ 2
    delta_lower_internal_less <-
      (1 - m_pre) * delta_lower / (0.25 + delta_upper ^ 2)
    delta_upper_internal_less <-
      (1 - m_pre) * delta_upper / delta_lower ^ 2
  } else {
    delta_lower_internal_greater <- delta_lower
    delta_upper_internal_greater <- delta_upper
    delta_lower_internal_less <- delta_lower
    delta_upper_internal_less <- delta_upper
  }
  
  # Load psi_fn list
  if (family == "Normal") {
    psi_fn_list <- generate_sub_G_fn(1)
    psi_fn_list_less <- generate_sub_G_fn(1)
  } else if (family == "Ber") {
    psi_fn_list <- generate_sub_B_fn(m_pre)
    psi_fn_list_less <- generate_sub_B_fn(1 - m_pre)
  } else if (family == "Bounded") {
    psi_fn_list <- generate_sub_E_fn()
    psi_fn_list_less <- generate_sub_E_fn()
  }
  
  # Compute weights and lambdas parameters
  if (alternative == "greater") {
    base_param <- compute_baseline(
      alpha,
      delta_lower_internal_greater,
      delta_upper_internal_greater,
      psi_fn_list,
      1,
      k_max
    )
    weights <- base_param$omega
    lambdas <- base_param$lambda
  } else if (alternative == "less") {
    base_param_less <- compute_baseline(
      alpha,
      delta_lower_internal_less,
      delta_upper_internal_less,
      psi_fn_list_less,
      1,
      k_max
    )
    weights <- base_param_less$omega
    if (family == "Normal" || family == "Ber") {
      lambdas <- -base_param_less$lambda
    } else if (family == "Bounded") {
      lambdas <- -m_pre * base_param_less$lambda / (1 - m_pre)
    }
  } else {
    base_param <- compute_baseline(
      alpha,
      delta_lower_internal_greater,
      delta_upper_internal_greater,
      psi_fn_list,
      1,
      k_max
    )
    base_param_less <- compute_baseline(
      alpha,
      delta_lower_internal_less,
      delta_upper_internal_less,
      psi_fn_list_less,
      1,
      k_max
    )
    weights <-
      c(base_param$omega / 2, base_param_less$omega / 2)
    if (family == "Normal" || family == "Ber") {
      lambdas <- c(base_param$lambda, -base_param_less$lambda)
    } else if (family == "Bounded") {
      lambdas <-
        c(base_param$lambda,
          -m_pre * base_param_less$lambda / (1 - m_pre))
    }
  }
  return(list(lambdas = lambdas, weights = weights))
}
