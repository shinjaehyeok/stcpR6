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

RCPP_MODULE(GLRCUNormalGreaterEx) {
  using namespace stcp;
  using GL = NormalGLRGreater;
  using GE = GLRCU<GL>;
  
  Rcpp::class_<Stcp<GE>>("GLRCUNormalGreaterBase")
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
  
  
  Rcpp::class_<GLRCUNormal<GL>>("GLRCUNormalGreater")
    .derives<Stcp<GE>>("GLRCUNormalGreaterBase")
    .constructor()
    .constructor<double, double, double, int>()
  ;
  
}

RCPP_MODULE(GLRCUNormalLessEx) {
  using namespace stcp;
  using GL = NormalGLRLess;
  using GE = GLRCU<GL>;
  
  Rcpp::class_<Stcp<GE>>("GLRCUNormalLessBase")
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
  
  
  Rcpp::class_<GLRCUNormal<GL>>("GLRCUNormalLess")
    .derives<Stcp<GE>>("GLRCUNormalLessBase")
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

RCPP_MODULE(GLRCUBerGreaterEx) {
  using namespace stcp;
  using GL = BerGLRGreater;
  using GE = GLRCU<GL>;
  
  Rcpp::class_<Stcp<GE>>("GLRCUBerGreaterBase")
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
  
  
  Rcpp::class_<GLRCUBer<GL>>("GLRCUBerGreater")
    .derives<Stcp<GE>>("GLRCUBerGreaterBase")
    .constructor()
    .constructor<double, double, int>()
  ;
  
}

RCPP_MODULE(GLRCUBerLessEx) {
  using namespace stcp;
  using GL = BerGLRLess;
  using GE = GLRCU<GL>;
  
  Rcpp::class_<Stcp<GE>>("GLRCUBerLessBase")
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
  
  
  Rcpp::class_<GLRCUBer<GL>>("GLRCUBerLess")
    .derives<Stcp<GE>>("GLRCUBerLessBase")
    .constructor()
    .constructor<double, double, int>()
  ;
  
}