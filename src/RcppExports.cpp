// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_StcpMixESTNormalEx();
RcppExport SEXP _rcpp_module_boot_StcpMixESRNormalEx();
RcppExport SEXP _rcpp_module_boot_StcpMixECUNormalEx();
RcppExport SEXP _rcpp_module_boot_StcpMixESTBerEx();
RcppExport SEXP _rcpp_module_boot_StcpMixESRBerEx();
RcppExport SEXP _rcpp_module_boot_StcpMixECUBerEx();
RcppExport SEXP _rcpp_module_boot_StcpMixESTBoundedrEx();
RcppExport SEXP _rcpp_module_boot_StcpMixESRBoundedrEx();
RcppExport SEXP _rcpp_module_boot_StcpMixECUBoundedrEx();
RcppExport SEXP _rcpp_module_boot_GLRCUNormalEx();
RcppExport SEXP _rcpp_module_boot_GLRCUBerEx();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_StcpMixESTNormalEx", (DL_FUNC) &_rcpp_module_boot_StcpMixESTNormalEx, 0},
    {"_rcpp_module_boot_StcpMixESRNormalEx", (DL_FUNC) &_rcpp_module_boot_StcpMixESRNormalEx, 0},
    {"_rcpp_module_boot_StcpMixECUNormalEx", (DL_FUNC) &_rcpp_module_boot_StcpMixECUNormalEx, 0},
    {"_rcpp_module_boot_StcpMixESTBerEx", (DL_FUNC) &_rcpp_module_boot_StcpMixESTBerEx, 0},
    {"_rcpp_module_boot_StcpMixESRBerEx", (DL_FUNC) &_rcpp_module_boot_StcpMixESRBerEx, 0},
    {"_rcpp_module_boot_StcpMixECUBerEx", (DL_FUNC) &_rcpp_module_boot_StcpMixECUBerEx, 0},
    {"_rcpp_module_boot_StcpMixESTBoundedrEx", (DL_FUNC) &_rcpp_module_boot_StcpMixESTBoundedrEx, 0},
    {"_rcpp_module_boot_StcpMixESRBoundedrEx", (DL_FUNC) &_rcpp_module_boot_StcpMixESRBoundedrEx, 0},
    {"_rcpp_module_boot_StcpMixECUBoundedrEx", (DL_FUNC) &_rcpp_module_boot_StcpMixECUBoundedrEx, 0},
    {"_rcpp_module_boot_GLRCUNormalEx", (DL_FUNC) &_rcpp_module_boot_GLRCUNormalEx, 0},
    {"_rcpp_module_boot_GLRCUBerEx", (DL_FUNC) &_rcpp_module_boot_GLRCUBerEx, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_stcpR6(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
