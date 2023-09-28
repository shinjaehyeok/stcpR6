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

# Modules from stcp_glr_export.cpp
Rcpp::loadModule(module = "GLRCUNormalEx", TRUE)
Rcpp::loadModule(module = "GLRCUBerEx", TRUE)

#' R6 Class representing a person
#'
#' A person has a name and a hair color.
#' @importFrom R6 R6Class
Stcp <- R6::R6Class("Stcp",
                      public = list(
                        #' @field name First or full name of the person.
                        name = NULL,
                        
                        #' @field hair Hair color of the person.
                        hair = NULL,
                        
                        #' @description
                        #' Create a new person object.
                        #' @param name Name.
                        #' @param hair Hair color.
                        #' @return A new `Person` object.
                        initialize = function(name = NA, hair = NA) {
                          self$name <- name
                          self$hair <- hair
                          self$greet()
                          private$stcpCpp <- NULL
                        },
                        
                        #' @description
                        #' Change hair color.
                        #' @param val New hair color.
                        set_hair = function(val) {
                          self$hair <- val
                        },
                        
                        #' @description
                        #' Say hi.
                        greet = function() {
                          cat(paste0("Hello, my name is ", self$name, ".\n"))
                        }
                      ),
                      private = list(
                        stcpCpp = NULL
                      )
)

# StcpNormal
# StcpBer
# StcpBounded


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
