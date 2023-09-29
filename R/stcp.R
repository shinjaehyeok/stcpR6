#' @title Stcp Class
#'
#' @description
#' Stcp class supports a unified framework for sequential tests and change
#' detection algorithms for streams of univariate (sub-)Gaussian, binary,
#'  and bounded random variables.
#'
#' @export
#' @importFrom R6 R6Class
Stcp <- R6::R6Class(
  "Stcp",
  public = list(
    #' @description
    #' Create a new Stcp object.
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
    #'
    #' @return A new `Person` object.
    initialize = function(method = c("ST", "SR", "CU", "GLRCU"),
                          family = c("Normal", "Ber", "Bounded"),
                          alternative = c("two.sided", "greater", "less"),
                          threshold = log(1 / 0.05),
                          m_pre = 0,
                          delta_lower = 0.1,
                          delta_upper = NULL,
                          k_max = 1000) {
      # Check input parameters
      method <- match.arg(method)
      family <- match.arg(family)
      alternative <- match.arg(alternative)
      
      if (threshold <= 0) {
        stop("threshold must be positive.")
      }
      
      alpha <- exp(-threshold)
      
      if (delta_lower <= 0 && method != "GLRCU") {
        stop("delta_lower must be positive.")
      }
      
      if (family == "Normal") {
        psi_fn_list <- generate_sub_G_fn(1)
      } else if (family == "Ber") {
        psi_fn_list <- generate_sub_B_fn(m_pre)
        if (alternative != "less") {
          if (m_pre + delta_lower > 1) {
            stop("The minimum of alternative / post-change parameter (m_pre + delta_lower) is greater than 1.")
          }
        }
        if (alternative != "greater") {
          if (m_pre - delta_lower < 0) {
            stop("The maximum of alternative / post-change parameter (m_pre - delta_lower) is less than 0.")
          }
        }
      } else if (family == "Bounded") {
        psi_fn_list <- generate_sub_E_fn()
        if (alternative == "two.sided") {
          # Note if alternative == "greater", bounded method technically does not require the upper bound 1.
          if (m_pre + delta_lower > 1) {
            stop("The minimum of alternative / post-change parameter (m_pre + delta_lower) is greater than 1.")
          }
        }
        if (alternative != "greater") {
          if (m_pre - delta_lower < 0) {
            stop("The maximum of alternative / post-change parameter (m_pre - delta_lower) is less than 0.")
          }
        }
      } else {
        stop("Unknown family")
      }
      
      if (is.null(delta_upper)) {
        # If delta_upper is NULL, we pick a reasonably large one
        delta_init_upper <-
          psi_fn_list$psi_star_inv(log(1 / alpha))
        
        if (delta_init_upper <= delta_lower) {
          delta_upper <- delta_lower
        } else {
          baseline_init <- compute_baseline(alpha,
                                            delta_lower,
                                            delta_init_upper,
                                            psi_fn_list,
                                            1,
                                            k_max)
          g_alpha <- baseline_init$g_alpha
          
          delta_upper <-
            max(delta_lower, psi_fn_list$psi_star_inv(g_alpha))
          
        }
        
      } else if (delta_upper <  delta_lower) {
        stop("If not NULL, delta_upper must be greater than or equal to delta_lower.")
      }
      
      if (family == "Bounded") {
        # Bounded family uses sub-E class internally 
        # So we convert it into the sub-E space.
        # where delta_E = m * delta / (sigma^2 + delta^2)
        delta_lower_internal <- m_pre * delta_lower / (0.25 + delta_upper^2)
        delta_upper_internal <- m_pre * delta_upper / delta_lower^2
      } else {
        delta_lower_internal <- delta_lower
        delta_upper_internal <- delta_upper
      }
      
      # Compute baseline parameters
      # TODO make two-sided and less version
      base_param <- compute_baseline(
        alpha,
        delta_lower_internal,
        delta_upper_internal,
        psi_fn_list,
        v_min,
        k_max
      )
      
      
      private$m_method <- method
      private$m_family <- family
      private$m_alternative <- alternative
      private$m_threshold <- threshold
      private$m_m_pre <- m_pre
      private$m_delta_lower <- delta_lower
      private$m_delta_upper <- delta_upper
      private$m_k_max <- k_max
      private$m_base_param <- base_param
    },
    greet = function() {
      cat(paste0("Hello, my name is ", self$name, ".\n"))
    }
  ),
  private = list(
    m_method = NULL,
    m_family = NULL,
    m_alternative = NULL,
    m_threshold = NULL,
    m_m_pre = NULL,
    m_delta_lower = NULL,
    m_delta_upper = NULL,
    m_k_max = NULL,
    m_base_param = NULL,
    m_stcpCpp = NULL
  )
)
