// stcp_export.cpp

// Include Rcpp system header file (e.g. <>)
#include <Rcpp.h>

// Include our definition of the student file (e.g. "")
#include "stcp.h"

// Expose (some of) the Student class
RCPP_MODULE(MixSTNormalEx) {
  using namespace stcp; // Name used to "loadModule" in R script
  using GE = ST<Normal>;
  Rcpp::class_<MixE<GE>>("MixSTNormal")       // This must be the C++ class name.
  .constructor()
  .constructor<std::vector<double>, std::vector<double>, std::vector<double>, double, double>()
  .method("getLogValue", &MixE<GE>::getLogValue)
  ;
}