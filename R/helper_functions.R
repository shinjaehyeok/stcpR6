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