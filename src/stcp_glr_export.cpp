// stcp_export.cpp

#include <Rcpp.h>

#include "stcp_export.h"

RCPP_MODULE(GLRCUNormalEx) {
  using namespace stcp;
  using GL = NormalGLR;
  using GE = GLRCU<GL>;
  
  Rcpp::class_<Stcp<GE>>("GLRCUNormalBase")
    .constructor()
  
  .method("getLogValue", &Stcp<GE>::getLogValue)
  .method("getThreshold", &Stcp<GE>::getThreshold)
  .method("isStopped", &Stcp<GE>::isStopped)
  .method("getTime", &Stcp<GE>::getTime)
  .method("getStoppedTime", &Stcp<GE>::getStoppedTime)
  .method("reset", &Stcp<GE>::reset)
  .method("updateLogValues", &Stcp<GE>::updateLogValues)
  .method("updateLogValuesUntilStop", &Stcp<GE>::updateLogValuesUntilStop)
  .method("updateAndReturnHistories", &Stcp<GE>::updateAndReturnHistories)
  ;
  
  
  Rcpp::class_<GLRCUNormal<GL>>("GLRCUNormal")
    .derives<Stcp<GE>>("GLRCUNormalBase")
    .constructor()
    .constructor<double, double, double, int>()
  ;
  
}

RCPP_MODULE(GLRCUBerEx) {
  using namespace stcp;
  using GL = BerGLR;
  using GE = GLRCU<GL>;
  
  Rcpp::class_<Stcp<GE>>("GLRCUBerBase")
    .constructor()
  
  .method("getLogValue", &Stcp<GE>::getLogValue)
  .method("getThreshold", &Stcp<GE>::getThreshold)
  .method("isStopped", &Stcp<GE>::isStopped)
  .method("getTime", &Stcp<GE>::getTime)
  .method("getStoppedTime", &Stcp<GE>::getStoppedTime)
  .method("reset", &Stcp<GE>::reset)
  .method("updateLogValues", &Stcp<GE>::updateLogValues)
  .method("updateLogValuesUntilStop", &Stcp<GE>::updateLogValuesUntilStop)
  .method("updateAndReturnHistories", &Stcp<GE>::updateAndReturnHistories)
  ;
  
  
  Rcpp::class_<GLRCUBer<GL>>("GLRCUBer")
    .derives<Stcp<GE>>("GLRCUBerBase")
    .constructor()
    .constructor<double, double, int>()
  ;
  
}