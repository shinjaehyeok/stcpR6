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
    .constructor<std::vector<double>,
                 std::vector<double>,
                 double,
                 double,
                 double>()
    ;
}
