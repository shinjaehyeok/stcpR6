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
    #' @param weights If not null, the input weights will be used to initialize Stcp object.
    #'
    #' @param lambdas If not null, the input lambdas will be used to initialize Stcp object.
    #'
    #' @param k_max Positive integer to determine the maximum number of baselines.
    #' For GLRCU method, it is used as the lookup window size for GLRCU statistics.
    #'
    #' @return A new `Person` object.
    initialize = function(method = c("ST", "SR", "CU", "GLRCU"),
                          family = c("Normal", "Ber", "Bounded"),
                          alternative = c("two.sided", "greater", "less"),
                          threshold = log(1 / 0.05),
                          m_pre = 0,
                          delta_lower = 0.1,
                          delta_upper = NULL,
                          weights = NULL,
                          lambdas = NULL,
                          k_max = 1000) {
      # Check input parameters
      method <- match.arg(method)
      family <- match.arg(family)
      alternative <- match.arg(alternative)
      
      if (threshold <= 0) {
        stop("threshold must be positive.")
      }
      
      alpha <- exp(-threshold)
      
      # If method = GLRCU, we do not use delta_lower and delta_upper
      if (method == "GLRCU") {
        if (family == "Normal") {
          if (alternative == "two.sided") {
            private$m_stcpCpp <- GLRCUNormal$new(threshold,
                                                 m_pre,
                                                 1,
                                                 k_max)
          } else if (alternative == "greater") {
            private$m_stcpCpp <- GLRCUNormalGreater$new(threshold,
                                                        m_pre,
                                                        1,
                                                        k_max)
          } else {
            private$m_stcpCpp <- GLRCUNormalLess$new(threshold,
                                                     m_pre,
                                                     1,
                                                     k_max)
          }
        } else if (family == "Ber") {
          if (alternative == "two.sided") {
            private$m_stcpCpp <- GLRCUBer$new(threshold,
                                              m_pre,
                                              k_max)
            
          } else if (alternative == "greater") {
            private$m_stcpCpp <- GLRCUBerGreater$new(threshold,
                                                     m_pre,
                                                     k_max)
          } else {
            private$m_stcpCpp <- GLRCUBerLess$new(threshold,
                                                  m_pre,
                                                  k_max)
          }
        } else {
          stop("Unsupported family for GLRCU method")
        }
        delta_lower <- 0
        weights <- NULL
        lambdas <- NULL
      } else {
        if (delta_lower <= 0 && method != "GLRCU") {
          stop("delta_lower must be positive.")
        }
        
        # Check delta_lower is within boundary for Ber and Bounded cases
        if (family == "Ber") {
          if (alternative != "less") {
            if (m_pre + delta_lower >= 1) {
              stop(
                "The minimum of alternative / post-change parameter (m_pre + delta_lower) is greater than or equal to 1."
              )
            }
          }
          if (alternative != "greater") {
            if (m_pre - delta_lower <= 0) {
              stop(
                "The maximum of alternative / post-change parameter (m_pre - delta_lower) is less than or equal to 0."
              )
            }
          }
        } else if (family == "Bounded") {
          if (alternative == "two.sided") {
            # Note if alternative == "greater", bounded method technically does not require the upper bound 1.
            if (m_pre + delta_lower > 1) {
              stop(
                "The minimum of alternative / post-change parameter (m_pre + delta_lower) is greater than 1."
              )
            }
          }
          if (alternative != "greater") {
            if (m_pre - delta_lower < 0) {
              stop(
                "The maximum of alternative / post-change parameter (m_pre - delta_lower) is less than 0."
              )
            }
          }
        }
        
        if (is.null(delta_upper)) {
          # If delta_upper is NULL, we pick a reasonably large one
          if (family == "Normal") {
            delta_upper_ <- 5
          } else {
            if (alternative == "greater") {
              delta_upper_ <- 1 - m_pre - 0.001
            } else if (alternative == "less") {
              delta_upper_ <- m_pre - 0.001
            } else {
              delta_upper_ <- min(1 - m_pre, m_pre) - 0.001
            }
          }
          
          delta_upper <- max(delta_lower, delta_upper_)
          
        } else if (delta_upper <  delta_lower) {
          stop("If not NULL, delta_upper must be greater than or equal to delta_lower.")
        }
        
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
        if (!is.null(weights) & !is.null(lambdas)) {
          if (length(weights) != length(lambdas)) {
            stop("Lengths of weights and lambdas are not same.")
          }
          if (length(weights) > k_max) {
            stop("Length of weights and lambdas exceed k_max.")
          }
        } else {
          if (alternative == "greater") {
            base_param <- compute_baseline(alpha,
                                           delta_lower_internal_greater,
                                           delta_upper_internal_greater,
                                           psi_fn_list,
                                           1,
                                           k_max)
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
            base_param <- compute_baseline(alpha,
                                           delta_lower_internal_greater,
                                           delta_upper_internal_greater,
                                           psi_fn_list,
                                           1,
                                           k_max)
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
              lambdas <- c(base_param$lambda,-base_param_less$lambda)
            } else if (family == "Bounded") {
              lambdas <-
                c(base_param$lambda,
                  -m_pre * base_param_less$lambda / (1 - m_pre))
            }
          }
          
        }
        
        
        # Initialize stcp Cpp-object
        if (family == "Normal") {
          if (method == "ST") {
            private$m_stcpCpp <- StcpMixESTNormal$new(threshold,
                                                      weights,
                                                      lambdas,
                                                      m_pre,
                                                      1)
          } else if (method == "SR") {
            private$m_stcpCpp <- StcpMixESRNormal$new(threshold,
                                                      weights,
                                                      lambdas,
                                                      m_pre,
                                                      1)
            
            
          } else if (method == "CU") {
            private$m_stcpCpp <- StcpMixECUNormal$new(threshold,
                                                      weights,
                                                      lambdas,
                                                      m_pre,
                                                      1)
            
          }
        } else if (family == "Ber") {
          if (method == "ST") {
            private$m_stcpCpp <- StcpMixESTBer$new(threshold,
                                                   weights,
                                                   lambdas,
                                                   m_pre)
          } else if (method == "SR") {
            private$m_stcpCpp <- StcpMixESRBer$new(threshold,
                                                   weights,
                                                   lambdas,
                                                   m_pre)
          } else if (method == "CU") {
            private$m_stcpCpp <- StcpMixECUBer$new(threshold,
                                                   weights,
                                                   lambdas,
                                                   m_pre)
          }
        } else if (family == "Bounded") {
          if (method == "ST") {
            private$m_stcpCpp <- StcpMixESTBounded$new(threshold,
                                                       weights,
                                                       lambdas,
                                                       m_pre)
          } else if (method == "SR") {
            private$m_stcpCpp <- StcpMixESRBounded$new(threshold,
                                                       weights,
                                                       lambdas,
                                                       m_pre)
          } else if (method == "CU") {
            private$m_stcpCpp <- StcpMixECUBounded$new(threshold,
                                                       weights,
                                                       lambdas,
                                                       m_pre)
          }
        }
      }
      private$m_method <- method
      private$m_family <- family
      private$m_alternative <- alternative
      private$m_threshold <- threshold
      private$m_alpha <- alpha
      private$m_m_pre <- m_pre
      private$m_delta_lower <- delta_lower
      private$m_delta_upper <- delta_upper
      private$m_k_max <- k_max
      private$m_weights <- weights
      private$m_lambdas <- lambdas
      
    },
    #' @description
    #' Print summary of Stcp object.
    print = function() {
      cat("stcp Model:\n")
      cat("- Method: ", private$m_method, "\n")
      cat("- Family: ", private$m_family, "\n")
      cat("- Alternative: ", private$m_alternative, "\n")
      cat("- Alpha: ", private$m_alpha, "\n")
      cat("- m_pre: ", private$m_m_pre, "\n")
      cat("- Num. of mixing components: ",
          length(private$m_weights),
          "\n")
      cat("- Obs. have been passed: ", self$getTime(), "\n")
      cat("- Current log value: ", self$getLogValue(), "\n")
      cat("- Is stopped before: ", self$isStopped(), "\n")
      cat("- Stopped time: ", self$getStoppedTime(), "\n")
    },
    getWeights = function() {
      private$m_weights
    },
    getLambdas = function() {
      private$m_lambdas
    },
    getLogValue = function() {
      private$m_stcpCpp$getLogValue()
    },
    getThreshold = function() {
      private$m_stcpCpp$getThreshold()
    },
    isStopped = function() {
      private$m_stcpCpp$isStopped()
    },
    getTime = function() {
      private$m_stcpCpp$getTime()
    },
    getStoppedTime = function() {
      private$m_stcpCpp$getStoppedTime()
    },
    reset = function() {
      private$m_stcpCpp$reset()
    },
    updateLogValues = function(xs) {
      private$m_stcpCpp$updateLogValues(xs)
    },
    updateLogValuesUntilStop = function(xs) {
      private$m_stcpCpp$updateLogValuesUntilStop(xs)
    },
    updateAndReturnHistories = function(xs) {
      private$m_stcpCpp$updateAndReturnHistories(xs)
    }
  ),
  private = list(
    m_method = NULL,
    m_family = NULL,
    m_alternative = NULL,
    m_threshold = NULL,
    m_alpha = NULL,
    m_m_pre = NULL,
    m_delta_lower = NULL,
    m_delta_upper = NULL,
    m_k_max = NULL,
    m_weights = NULL,
    m_lambdas = NULL,
    m_stcpCpp = NULL
  ),
  cloneable = FALSE
)
