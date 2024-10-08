#' @title NormalCS Class
#'
#' @description
#' NormalCS class is used to compute always-valid confidence sequence
#' for the standard normal process based on the STCP method.
#'
#' @export
#' @importFrom R6 R6Class
#' 
#' @examples
#' # Initialize two-sided standard normal confidence sequence
#' # optimized for the interval [10, 100]
#' normal_cs <- NormalCS$new(
#'   alternative = "two.sided",
#'   alpha = 0.05,
#'   n_upper = 100,
#'   n_lower = 10
#'   )
#'   
#' # Compute confidence interval at n = 20 when observed sample mean = 0.5
#' normal_cs$computeInterval(20, x_bar = 0.5)
#' 
#' # (Advanced) NormalCS supports general variance process.
#' # Both n_upper and n_lower can be general positive numbers.
#' normal_cs2 <- NormalCS$new(
#'   alternative = "two.sided",
#'   alpha = 0.05,
#'   n_upper = 100.5,
#'   n_lower = 10.5
#'   )
#' # Confidence interval at n = 20.5
#' normal_cs$computeInterval(20.5, x_bar = 0.5)
#' 
NormalCS <- R6::R6Class(
  "NormalCS",
  public = list(
    #' @description
    #' Create a new NormalCS object.
    #'
    #' @param alternative Alternative / post-change mean space
    #' * two.sided: Two-sided test / change detection
    #' * greater: Alternative /post-change mean is greater than null / pre-change one
    #' * less:  Alternative /post-change mean is less than null / pre-change one
    #'
    #' @param alpha Upper bound on the type 1 error of the confidence sequence.
    #'
    #' @param n_upper Upper bound of the target sample interval
    #' 
    #' @param n_lower Lower bound of the target sample interval
    #'
    #' @param weights If not null, the input weights will be used to initialize the object
    #' instead of \code{n_upper} and \code{n_lower}.
    #'
    #' @param lambdas If not null, the input lambdas will be used to initialize the object.
    #' instead of \code{n_upper} and \code{n_lower}.
    #' 
    #' @param skip_g_alpha If true, we do not compute g_alpha and use log(1/alpha) instead.
    #' 
    #' @param k_max Positive integer to determine the maximum number of baselines.
    #'
    #' @return A new `NormalCS` object.
    #' 
    initialize = function(alternative = c("two.sided", "greater", "less"),
                          alpha = 0.05,
                          n_upper = 1000,
                          n_lower = 1,
                          weights = NULL,
                          lambdas = NULL,
                          skip_g_alpha = TRUE,
                          k_max = 1000) {
      # Check input parameters
      alternative <- match.arg(alternative)
      
      if (alpha <= 0 | alpha > 1) {
        stop("alpha must be strictly inbetween 0 and 1.")
      }
       
      # Compute weights and lambdas parameters
      if (!is.null(weights) & !is.null(lambdas)) {
        if (length(weights) != length(lambdas)) {
          stop("Lengths of weights and lambdas are not same.")
        }
        if (length(weights) > k_max) {
          stop("Length of weights and lambdas exceed k_max.")
        }
        stcp <- Stcp$new(
          method = "ST",
          family = "Normal",
          alternative = alternative,
          threshold = log(1 / alpha),
          m_pre = 0,
          weights = weights,
          lambdas = lambdas,
          k_max = k_max
        )
      } else {
        base_param <- compute_baseline_for_sample_size(alpha,
                                                       n_upper,
                                                       n_lower,
                                                       generate_sub_G_fn(),
                                                       skip_g_alpha,
                                                       n_lower,
                                                       k_max)
        stcp <- Stcp$new(
          method = "ST",
          family = "Normal",
          alternative = alternative,
          threshold = log(1 / alpha),
          m_pre = 0,
          delta_lower = base_param$delta_lower,
          delta_upper = base_param$delta_upper,
          k_max = k_max
        )
        weights <- stcp$getWeights()
        lambdas <- stcp$getLambdas()
      }
      
      private$m_stcp <- stcp
      private$m_alternative <- alternative
      private$m_alpha <- alpha
      private$m_n_upper <- n_upper
      private$m_n_lower <- n_lower
      private$m_k_max <- k_max
      private$m_weights <- weights
      private$m_lambdas <- lambdas
    },
    #' @description
    #' Print summary of Stcp object.
    print = function() {
      cat("NormalCS Object:\n")
      cat("- Alternative: ", private$m_alternative, "\n")
      cat("- Alpha: ", private$m_alpha, "\n")
      cat("- Num. of mixing components: ",
          length(private$m_weights),
          "\n")
    },
    #' @description
    #' Return the upper bound on the type 1 error
    getAlpha = function() {
      private$m_alpha
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
    #' Compute the width of confidence interval at time n.
    #' 
    #' @param n Positive time.
    computeWidth = function(n) {
      compute_gap_from_boundary <- function(mu){
        private$m_stcp$reset()
        # log value at x_bar = 0, n and mu.
        log_value <- private$m_stcp$updateAndReturnHistoriesByAvgs(-mu, n)
        return(log_value - private$m_stcp$getThreshold())
      }
      chernoff_width <- sqrt(2 * log(1/private$m_alpha) / n)
      if (private$m_alternative == "greater") {
        search_l <- -20 * chernoff_width
        search_u <- -chernoff_width
      } else {
        search_l <- chernoff_width
        search_u <- 20 * chernoff_width
      }
      boundary_search <- stats::uniroot(compute_gap_from_boundary,
                                        c(search_l, search_u), tol = 1e-10)
      return(abs(boundary_search$root))
    },
    #' @description
    #' Compute a vector of two end points of confidence interval
    #' at time n 
    #' 
    #' @param n Positive time.
    #' @param x_bar The center of the confidence interval.
    computeInterval = function(n, x_bar = 0) {
      width <- self$computeWidth(n)
      if (private$m_alternative == "greater") {
        interval_out <- c(x_bar - width, Inf)
      } else if (private$m_alternative == "less") {
        interval_out <- c(-Inf, x_bar + width)
      } else {
        interval_out <- c(x_bar - width, x_bar + width)
      }
      return(interval_out)
    }
  ),
  private = list(
    m_alternative = NULL,
    m_alpha = NULL,
    m_n_upper = NULL,
    m_n_lower = NULL,
    m_k_max = NULL,
    m_weights = NULL,
    m_lambdas = NULL,
    m_stcp = NULL
  ),
  cloneable = FALSE
)
