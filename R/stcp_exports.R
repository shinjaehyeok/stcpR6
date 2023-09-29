# Modules from stcp_export.cpp
Rcpp::loadModule(module = "StcpMixESTNormalEx", TRUE)
Rcpp::loadModule(module = "StcpMixESRNormalEx", TRUE)
Rcpp::loadModule(module = "StcpMixECUNormalEx", TRUE)

Rcpp::loadModule(module = "StcpMixESTBerEx", TRUE)
Rcpp::loadModule(module = "StcpMixESRBerEx", TRUE)
Rcpp::loadModule(module = "StcpMixECUBerEx", TRUE)

Rcpp::loadModule(module = "StcpMixESTBoundedrEx", TRUE)
Rcpp::loadModule(module = "StcpMixESRBoundedrEx", TRUE)
Rcpp::loadModule(module = "StcpMixECUBoundedrEx", TRUE)

# Modules from stcp_glrcu_export.cpp
Rcpp::loadModule(module = "GLRCUNormalEx", TRUE)
Rcpp::loadModule(module = "GLRCUNormalGreaterEx", TRUE)
Rcpp::loadModule(module = "GLRCUNormalLessEx", TRUE)

Rcpp::loadModule(module = "GLRCUBerEx", TRUE)
Rcpp::loadModule(module = "GLRCUBerGreaterEx", TRUE)
Rcpp::loadModule(module = "GLRCUBerLessEx", TRUE)

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




#' Initialize an object of GLRCUBer class
#'
#' @param threshold stopping threshold
#' @param p mean of H0
#' @param window_size Window size for GLR computation.
#'
#' @return An object of StcpGLRCUBer class
#' @export
#'
#' @examples
makeGLRCUBer_ <- function(threshold = log(1 / 0.05),
                          p = 0,
                          window_size = 100L) {
  return(GLRCUBer$new(threshold,
                      p,
                      window_size))
}
