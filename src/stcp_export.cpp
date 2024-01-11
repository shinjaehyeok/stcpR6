// stcp_export.cpp

#include <Rcpp.h>

#include "stcp_export.h"

RCPP_MODULE(HelperEx) {
  using namespace stcp;
  Rcpp::function("logSumExp", &logSumExp,
           "Compute log-sum-exp of a numeric vector.");
}

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
    .method("updateLogValuesByAvgs", &Stcp<MixE<GE>>::updateLogValuesByAvgs)
    .method("updateLogValuesUntilStopByAvgs", &Stcp<MixE<GE>>::updateLogValuesUntilStopByAvgs)
    .method("updateAndReturnHistoriesByAvgs", &Stcp<MixE<GE>>::updateAndReturnHistoriesByAvgs)
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

RCPP_MODULE(StcpMixESRNormalEx) {
  using namespace stcp;
  using GE = SR<Normal>;
  
  Rcpp::class_<Stcp<MixE<GE>>>("StcpMixESRNormalBase")
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
  .method("updateLogValuesByAvgs", &Stcp<MixE<GE>>::updateLogValuesByAvgs)
  .method("updateLogValuesUntilStopByAvgs", &Stcp<MixE<GE>>::updateLogValuesUntilStopByAvgs)
  .method("updateAndReturnHistoriesByAvgs", &Stcp<MixE<GE>>::updateAndReturnHistoriesByAvgs)
  ;
  
  
  Rcpp::class_<StcpNormal<GE>>("StcpMixESRNormal")
    .derives<Stcp<MixE<GE>>>("StcpMixESRNormalBase")
    .constructor()
    .constructor<double, 
  std::vector<double>,
  std::vector<double>,
  double,
  double>()
    ;
}

RCPP_MODULE(StcpMixECUNormalEx) {
  using namespace stcp;
  using GE = CU<Normal>;
  
  Rcpp::class_<Stcp<MixE<GE>>>("StcpMixECUNormalBase")
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
  .method("updateLogValuesByAvgs", &Stcp<MixE<GE>>::updateLogValuesByAvgs)
  .method("updateLogValuesUntilStopByAvgs", &Stcp<MixE<GE>>::updateLogValuesUntilStopByAvgs)
  .method("updateAndReturnHistoriesByAvgs", &Stcp<MixE<GE>>::updateAndReturnHistoriesByAvgs)
  ;
  
  
  Rcpp::class_<StcpNormal<GE>>("StcpMixECUNormal")
    .derives<Stcp<MixE<GE>>>("StcpMixECUNormalBase")
    .constructor()
    .constructor<double, 
  std::vector<double>,
  std::vector<double>,
  double,
  double>()
    ;
}

RCPP_MODULE(StcpMixESTBerEx) {
  using namespace stcp;
  using GE = ST<Ber>;
  
  Rcpp::class_<Stcp<MixE<GE>>>("StcpMixESTBerBase")
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
  .method("updateLogValuesByAvgs", &Stcp<MixE<GE>>::updateLogValuesByAvgs)
  .method("updateLogValuesUntilStopByAvgs", &Stcp<MixE<GE>>::updateLogValuesUntilStopByAvgs)
  .method("updateAndReturnHistoriesByAvgs", &Stcp<MixE<GE>>::updateAndReturnHistoriesByAvgs)
  ;
  
  
  Rcpp::class_<StcpBer<GE>>("StcpMixESTBer")
    .derives<Stcp<MixE<GE>>>("StcpMixESTBerBase")
    .constructor()
    .constructor<double, 
  std::vector<double>,
  std::vector<double>,
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
  .method("updateLogValuesByAvgs", &Stcp<MixE<GE>>::updateLogValuesByAvgs)
  .method("updateLogValuesUntilStopByAvgs", &Stcp<MixE<GE>>::updateLogValuesUntilStopByAvgs)
  .method("updateAndReturnHistoriesByAvgs", &Stcp<MixE<GE>>::updateAndReturnHistoriesByAvgs)
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

RCPP_MODULE(StcpMixECUBerEx) {
  using namespace stcp;
  using GE = CU<Ber>;
  
  Rcpp::class_<Stcp<MixE<GE>>>("StcpMixECUBerBase")
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
  .method("updateLogValuesByAvgs", &Stcp<MixE<GE>>::updateLogValuesByAvgs)
  .method("updateLogValuesUntilStopByAvgs", &Stcp<MixE<GE>>::updateLogValuesUntilStopByAvgs)
  .method("updateAndReturnHistoriesByAvgs", &Stcp<MixE<GE>>::updateAndReturnHistoriesByAvgs)
  ;
  
  
  Rcpp::class_<StcpBer<GE>>("StcpMixECUBer")
    .derives<Stcp<MixE<GE>>>("StcpMixECUBerBase")
    .constructor()
    .constructor<double, 
  std::vector<double>,
  std::vector<double>,
  double>()
    ;
  
}

RCPP_MODULE(StcpMixESTBoundedrEx) {
  using namespace stcp;
  using GE = ST<Bounded>;
  
  Rcpp::class_<Stcp<MixE<GE>>>("StcpMixESTBoundedBase")
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
  
  
  Rcpp::class_<StcpBounded<GE>>("StcpMixESTBounded")
    .derives<Stcp<MixE<GE>>>("StcpMixESTBoundedBase")
    .constructor()
    .constructor<double, 
  std::vector<double>,
  std::vector<double>,
  double>()
    ;
  
}

RCPP_MODULE(StcpMixESRBoundedrEx) {
  using namespace stcp;
  using GE = SR<Bounded>;
  
  Rcpp::class_<Stcp<MixE<GE>>>("StcpMixESRBoundedBase")
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
  
  
  Rcpp::class_<StcpBounded<GE>>("StcpMixESRBounded")
    .derives<Stcp<MixE<GE>>>("StcpMixESRBoundedBase")
    .constructor()
    .constructor<double, 
                 std::vector<double>,
                 std::vector<double>,
                 double>()
    ;
  
}

RCPP_MODULE(StcpMixECUBoundedrEx) {
  using namespace stcp;
  using GE = CU<Bounded>;
  
  Rcpp::class_<Stcp<MixE<GE>>>("StcpMixECUBoundedBase")
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
  
  
  Rcpp::class_<StcpBounded<GE>>("StcpMixECUBounded")
    .derives<Stcp<MixE<GE>>>("StcpMixECUBoundedBase")
    .constructor()
    .constructor<double, 
  std::vector<double>,
  std::vector<double>,
  double>()
    ;
  
}
