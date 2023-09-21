#' Initialize an object of MixSTNormal class
#'
#' @param weights Vector of weights
#' @param log_values Vector of log-values
#' @param lambdas Vector of lambdas
#' @param mu0 mean of H0
#' @param sig sigma
#'  
#' @return An object of MixSTNormal class
#' @export
#'
#' @examples
makeMixSTNormal <- function(weights = 1,
                            log_values = 0,
                            lambdas = 0,
                            mu0 = 0,
                            sig = 1) {
  
  return(MixSTNormal$new(weights,
                         log_values,
                         lambdas,
                         mu0,
                         sig))
}