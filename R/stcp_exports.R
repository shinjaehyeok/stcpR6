Rcpp::loadModule(module = "StcpMixESTNormalEx", TRUE)

#' Initialize an object of StcpMixSTNormal class
#'
#' @param threshold stopping threshold
#' @param weights Vector of weights
#' @param lambdas Vector of lambdas
#' @param mu0 mean of H0
#' @param sig sigma
#'
#' @return An object of StcpMixSTNormal class
#' @export
#'
#' @examples
makeStcpMixESTNormal_ <- function(threshold = log(1 / 0.05),
                                  weights = 1,
                                  lambdas = 1,
                                  mu0 = 0,
                                  sig = 1) {
  return(StcpMixESTNormal$new(threshold, weights,
                              lambdas,
                              mu0,
                              sig))
}

Rcpp::loadModule(module = "StcpMixESRBerEx", TRUE)

#' Initialize an object of StcpMixSRBer class
#'
#' @param threshold stopping threshold
#' @param weights Vector of weights
#' @param lambdas Vector of lambdas
#' @param p mean of H0
#'
#' @return An object of StcpMixSRBer class
#' @export
#'
#' @examples
makeStcpMixESRBer_ <- function(threshold = log(1 / 0.05),
                               weights = 1,
                               lambdas = 1,
                               p = 0) {
  return(StcpMixESRBer$new(threshold, weights,
                           lambdas,
                           p))
}


Rcpp::loadModule(module = "StcpGLRCUBerEx", TRUE)

#' Initialize an object of StcpGLRCUBer class
#'
#' @param threshold stopping threshold
#' @param p mean of H0
#' @param window_size Window size for GLR computation.
#'
#' @return An object of StcpGLRCUBer class
#' @export
#'
#' @examples
makeStcpGLRCUBer_ <- function(threshold = log(1 / 0.05),
                              p = 0,
                              window_size = 100L) {
  return(StcpGLRCUBer$new(threshold,
                          p,
                          window_size))
}
