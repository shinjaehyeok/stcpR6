// stcp_export.cpp

#include <Rcpp.h>

#include "stcp_export.h"

RCPP_MODULE(StcpMixESTNormalEx) {
  using namespace stcp;
  using GE = ST<Normal>;

  Rcpp::class_<Stcp<MixE<GE>>>("StcpMixESTNormalBase")
    .constructor()
  
    .method("getLogValue", &Stcp<MixE<GE>>::getLogValue)
    .method("getThreshold", &Stcp<MixE<GE>>::getThreshold)
    .method("isStopped", &Stcp<MixE<GE>>::isStopped)
    .method("getTime", &Stcp<MixE<GE>>::getTime)
    .method("getStoppedTime", &Stcp<MixE<GE>>::getStoppedTime)
    .method("reset", &Stcp<MixE<GE>>::reset)
    .method("updateLogValues", &Stcp<MixE<GE>>::updateLogValues)
    .method("updateLogValuesUntilStop", &Stcp<MixE<GE>>::updateLogValuesUntilStop)
    .method("updateAndReturnHistories", &Stcp<MixE<GE>>::updateAndReturnHistories)
    ;

    
  Rcpp::class_<StcpNormal<GE>>("StcpMixESTNormal")
    .derives<Stcp<MixE<GE>>>("StcpMixESTNormalBase")
    .constructor()
    .constructor<double, 
                 std::vector<double>,
                 std::vector<double>,
                 double,
                 double>()
    ;
}

RCPP_MODULE(StcpMixESRBerEx) {
  using namespace stcp;
  using GE = SR<Ber>;
  
  Rcpp::class_<Stcp<MixE<GE>>>("StcpMixESRBerBase")
    .constructor()
  
  .method("getLogValue", &Stcp<MixE<GE>>::getLogValue)
  .method("getThreshold", &Stcp<MixE<GE>>::getThreshold)
  .method("isStopped", &Stcp<MixE<GE>>::isStopped)
  .method("getTime", &Stcp<MixE<GE>>::getTime)
  .method("getStoppedTime", &Stcp<MixE<GE>>::getStoppedTime)
  .method("reset", &Stcp<MixE<GE>>::reset)
  .method("updateLogValues", &Stcp<MixE<GE>>::updateLogValues)
  .method("updateLogValuesUntilStop", &Stcp<MixE<GE>>::updateLogValuesUntilStop)
  .method("updateAndReturnHistories", &Stcp<MixE<GE>>::updateAndReturnHistories)
  ;
  
  
  Rcpp::class_<StcpBer<GE>>("StcpMixESRBer")
    .derives<Stcp<MixE<GE>>>("StcpMixESRBerBase")
    .constructor()
    .constructor<double, 
                 std::vector<double>,
                 std::vector<double>,
                 double>()
    ;
  
}

RCPP_MODULE(StcpGLRCUBerEx) {
  using namespace stcp;
  using GE = GLRCU<BerGLR>;
  
  Rcpp::class_<Stcp<GE>>("StcpGLRCUBerBase")
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
  
  
  Rcpp::class_<StcpBer<GE>>("StcpGLRCUBer")
    .derives<Stcp<GE>>("StcpGLRCUBerBase")
    .constructor()
    .constructor<double, double, int>()
    ;
  
}
