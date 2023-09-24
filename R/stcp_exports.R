Rcpp::loadModule(module = "StcpMixESTNormalEx", TRUE)

#' Initialize an object of StcpMixSTNormal class
#'
#' @param weights Vector of weights
#' @param lambdas Vector of lambdas
#' @param mu0 mean of H0
#' @param sig sigma
#' @param threshold stopping threshold
#'
#' @return An object of StcpMixSTNormal class
#' @export
#'
#' @examples
makeStcpMixESTNormal_ <- function(weights = 1,
                                  lambdas = 1,
                                  mu0 = 0,
                                  sig = 1,
                                  threshold = log(1 / 0.05)) {
  
  return(StcpMixESTNormal$new(weights,
                              lambdas,
                              mu0,
                              sig,
                              threshold))
}
