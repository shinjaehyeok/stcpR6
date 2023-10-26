# Modules from stcp_export.cpp
Rcpp::loadModule(module = "HelperEx", TRUE)

Rcpp::loadModule(module = "StcpMixESTNormalEx", TRUE)
Rcpp::loadModule(module = "StcpMixESRNormalEx", TRUE)
Rcpp::loadModule(module = "StcpMixECUNormalEx", TRUE)

Rcpp::loadModule(module = "StcpMixESTBerEx", TRUE)
Rcpp::loadModule(module = "StcpMixESRBerEx", TRUE)
Rcpp::loadModule(module = "StcpMixECUBerEx", TRUE)

Rcpp::loadModule(module = "StcpMixESTBoundedrEx", TRUE)
Rcpp::loadModule(module = "StcpMixESRBoundedrEx", TRUE)
Rcpp::loadModule(module = "StcpMixECUBoundedrEx", TRUE)

# Modules from stcp_glrcu_export.cpp
Rcpp::loadModule(module = "GLRCUNormalEx", TRUE)
Rcpp::loadModule(module = "GLRCUNormalGreaterEx", TRUE)
Rcpp::loadModule(module = "GLRCUNormalLessEx", TRUE)

Rcpp::loadModule(module = "GLRCUBerEx", TRUE)
Rcpp::loadModule(module = "GLRCUBerGreaterEx", TRUE)
Rcpp::loadModule(module = "GLRCUBerLessEx", TRUE)
