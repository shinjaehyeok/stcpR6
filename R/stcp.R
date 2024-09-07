#' @title Stcp Class
#'
#' @description
#' Stcp class supports a unified framework for sequential tests and change
#' detection algorithms for streams of univariate (sub-)Gaussian, binary,
#'  and bounded random variables.
#'
#' @export
#' @importFrom R6 R6Class
#'
#' @examples
#' # Sequential Normal mean test H0: mu <= 0
#' # Initialize stcp object for this test.
#' stcp <- Stcp$new(method = "ST",
#'                  family = "Normal",
#'                  alternative = "greater",
#'                  threshold = log(1 / 0.05),
#'                  m_pre = 0)
#'
#' # Update the observations
#' obs <- c(1.0, 3.0, 2.0)
#' stcp$updateLogValuesUntilStop(obs)
#'
#' # Check whether the sequential test is stopped
#' stcp$isStopped() # TRUE
#'
#' # Check when the test was stopped
#' stcp$getStoppedTime() # 3
#'
#' # Although the number of obervaions was 4, the test was stopped at 3.
#' stcp$getTime() # 3
#'
#' # Get the log value of the mixutre of e-values at the current time (3)
#' stcp$getLogValue() # 4.425555
#'
#' # ...which is higher than the threshold log(1 / 0.05) ~ 2.996
#' stcp$getThreshold() # 2.995732
#'
#' # Reset the test object
#' stcp$reset()
#'
#' # Rerun the test but, at this time, we track updated log values
#' log_values <- stcp$updateAndReturnHistories(obs)
#' print(log_values) # 0.1159777 2.7002207 4.4255551 1.9746508
#'
#' # Again, the test was stopped at 3rd observation
#' stcp$getStoppedTime() # 3
#'
#' # But, at this time, log values were evaluated until the 4th observation.
#' stcp$getTime() # 4
#'
#' # Print overall summary
#' stcp # or stcp$print() or print(stcp)
#' # stcp Model:
#' #   - Method:  ST
#' # - Family:  Normal
#' # - Alternative:  greater
#' # - Alpha:  0.05
#' # - m_pre:  0
#' # - Num. of mixing components:  55
#' # - Obs. have been passed:  4
#' # - Current log value:  1.974651
#' # - Is stopped before:  TRUE
#' # - Stopped time:  3
#'
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
    #' @return A new `Stcp` object.
    #'
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
      
      delta_check_output <- checkDeltaRange(method, family, alternative, m_pre, delta_lower, delta_upper)
      
      if (delta_check_output$is_acceptable) {
        delta_upper <- delta_check_output$delta_upper
      } else {
        stop(delta_check_output$error_message)
      }
      
      # Initialize stcp for GLRCU method
      # If method = GLRCU, we do not use delta_lower and delta_upper
      if (method == "GLRCU") {
        if (family == "Normal") {
          if (alternative == "two.sided") {
            private$m_stcpCpp <- GLRCUNormal$new(threshold, m_pre, 1, k_max)
          } else if (alternative == "greater") {
            private$m_stcpCpp <- GLRCUNormalGreater$new(threshold, m_pre, 1, k_max)
          } else {
            private$m_stcpCpp <- GLRCUNormalLess$new(threshold, m_pre, 1, k_max)
          }
        } else if (family == "Ber") {
          if (alternative == "two.sided") {
            private$m_stcpCpp <- GLRCUBer$new(threshold, m_pre, k_max)
          } else if (alternative == "greater") {
            private$m_stcpCpp <- GLRCUBerGreater$new(threshold, m_pre, k_max)
          } else {
            private$m_stcpCpp <- GLRCUBerLess$new(threshold, m_pre, k_max)
          }
        } else {
          stop("Unsupported family for GLRCU method")
        }
        delta_lower <- 0
        weights <- NULL
        lambdas <- NULL
        
        return()
      }
      
      # Initialize stcp for non-GLR methods
      # If input weights or lamdas is empty then
      # non-GRL methods are initialized by lambdas and weights
      
      if (!is.null(lambdas)) {
        if (length(lambdas) > k_max) {
          stop("Length of lambdas exceed k_max.")
        }
        if (is.null(weights)) {
          # If user input lambdas but not specified weights then
          # we use the uniform weight by default
          weights <- rep(1.0 / length(lambdas), length(lambdas))
        } else {
          if (length(weights) != length(lambdas)) {
            stop("Lengths of weights and lambdas are not same.")
          }
        }
      } else {
        exp_params <- convertDeltaToExpParams(family,
                                              alternative,
                                              threshold,
                                              m_pre,
                                              delta_lower,
                                              delta_upper,
                                              k_max)
        weights <- exp_params$weights
        lambdas <- exp_params$lambdas
      }
      
      # Initialize stcp Cpp-object
      if (family == "Normal") {
        if (method == "ST") {
          private$m_stcpCpp <- StcpMixESTNormal$new(threshold, weights, lambdas, m_pre, 1)
        } else if (method == "SR") {
          private$m_stcpCpp <- StcpMixESRNormal$new(threshold, weights, lambdas, m_pre, 1)
          
          
        } else if (method == "CU") {
          private$m_stcpCpp <- StcpMixECUNormal$new(threshold, weights, lambdas, m_pre, 1)
          
        }
      } else if (family == "Ber") {
        if (method == "ST") {
          private$m_stcpCpp <- StcpMixESTBer$new(threshold, weights, lambdas, m_pre)
        } else if (method == "SR") {
          private$m_stcpCpp <- StcpMixESRBer$new(threshold, weights, lambdas, m_pre)
        } else if (method == "CU") {
          private$m_stcpCpp <- StcpMixECUBer$new(threshold, weights, lambdas, m_pre)
        }
      } else if (family == "Bounded") {
        if (method == "ST") {
          private$m_stcpCpp <- StcpMixESTBounded$new(threshold, weights, lambdas, m_pre)
        } else if (method == "SR") {
          private$m_stcpCpp <- StcpMixESRBounded$new(threshold, weights, lambdas, m_pre)
        } else if (method == "CU") {
          private$m_stcpCpp <- StcpMixECUBounded$new(threshold, weights, lambdas, m_pre)
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
    #' @description
    #' Return weights of mixture of e-values / e-detectors.
    getWeights = function() {
      private$m_weights
    },
    #' @description
    #' Return lambda parameters of mixture of e-values / e-detectors.
    getLambdas = function() {
      private$m_lambdas
    },
    #' @description
    #' Return the log value of mixture of e-values / e-detectors.
    getLogValue = function() {
      private$m_stcpCpp$getLogValue()
    },
    #' @description
    #' Return the threshold of the sequential test / change detection
    getThreshold = function() {
      private$m_stcpCpp$getThreshold()
    },
    #' @description
    #' Return TRUE if the sequential test / change detection was stopped by crossing the threshold.
    isStopped = function() {
      private$m_stcpCpp$isStopped()
    },
    #' @description
    #' Return the number of observations having been passed.
    getTime = function() {
      private$m_stcpCpp$getTime()
    },
    #' @description
    #' Return the stopped time. If it has been never stopped, return zero.
    getStoppedTime = function() {
      private$m_stcpCpp$getStoppedTime()
    },
    #' @description
    #' Reset the stcp object to the initial setup.
    reset = function() {
      private$m_stcpCpp$reset()
    },
    #' @description
    #' Update the log value and related fields by passing a vector of observations.
    #'
    #' @param xs A numeric vector of observations.
    updateLogValues = function(xs) {
      private$m_stcpCpp$updateLogValues(xs)
    },
    #' @description
    #' Update the log value and related fields until the log value is crossing the boundary.
    #'
    #' @param xs A numeric vector of observations.
    updateLogValuesUntilStop = function(xs) {
      private$m_stcpCpp$updateLogValuesUntilStop(xs)
    },
    #' @description
    #' Update the log value and related fields then return updated log values by passing a vector of observations.
    #'
    #' @param xs A numeric vector of observations.
    updateAndReturnHistories = function(xs) {
      private$m_stcpCpp$updateAndReturnHistories(xs)
    },
    #' @description
    #' Update the log value and related fields by passing
    #' a vector of averages and number of corresponding samples.
    #'
    #' @param x_bars A numeric vector of averages.
    #' @param ns A numeric vector of sample sizes.
    updateLogValuesByAvgs = function(x_bars, ns) {
      private$m_stcpCpp$updateLogValuesByAvgs(x_bars, ns)
    },
    #' @description
    #' Update the log value and related fields by passing
    #' a vector of averages and number of corresponding samples
    #' until the log value is crossing the boundary.
    #'
    #' @param x_bars A numeric vector of averages.
    #' @param ns A numeric vector of sample sizes.
    updateLogValuesUntilStopByAvgs = function(x_bars, ns) {
      private$m_stcpCpp$updateLogValuesUntilStopByAvgs(x_bars, ns)
    },
    #' @description
    #' Update the log value and related fields then return updated log values
    #' a vector of averages and number of corresponding samples.
    #'
    #' @param x_bars A numeric vector of averages.
    #' @param ns A numeric vector of sample sizes.
    updateAndReturnHistoriesByAvgs = function(x_bars, ns) {
      private$m_stcpCpp$updateAndReturnHistoriesByAvgs(x_bars, ns)
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
