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