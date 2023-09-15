// stcp_export.cpp

// Include Rcpp system header file (e.g. <>)
#include <Rcpp.h>

// Include our definition of the student file (e.g. "")
#include "stcp.h"

// Expose (some of) the Student class
RCPP_MODULE(MixEEx) {
  using namespace stcp; // Name used to "loadModule" in R script
  Rcpp::class_<MixE>("MixE")       // This must be the C++ class name.
  .constructor()
  .constructor<std::vector<double>>()
  .constructor<std::vector<double>, std::vector<double>>()
  .method("print", &MixE::print)
  .method("getLogMixedValue", &MixE::getLogMixedValue)
  ;
}