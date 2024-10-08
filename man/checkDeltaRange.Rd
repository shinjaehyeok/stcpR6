% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{checkDeltaRange}
\alias{checkDeltaRange}
\title{Check whether the input delta parameters are acceptable}
\usage{
checkDeltaRange(method, family, alternative, m_pre, delta_lower, delta_upper)
}
\arguments{
\item{method}{Method of the sequential procedure.
\itemize{
\item ST: Sequential test based on a mixture of E-values.
\item SR: Sequential change detection based on e-SR procedure.
\item CU: Sequential change detection based on e-CUSUM procedure.
\item GLRCU: Sequential change detection based on GLR-CUSUM procedure.
}}

\item{family}{Distribution of underlying univariate observations.
\itemize{
\item Normal: (sub-)Gaussian with sigma = 1.
\item Ber: Bernoulli distribution on \{0,1\}.
\item Bounded: General bounded distribution on [0,1]
}}

\item{alternative}{Alternative / post-change mean space
\itemize{
\item two.sided: Two-sided test / change detection
\item greater: Alternative /post-change mean is greater than null / pre-change one
\item less:  Alternative /post-change mean is less than null / pre-change one
}}

\item{m_pre}{The boundary of mean parameter in null / pre-change space}

\item{delta_lower}{Minimum gap between null / pre-change space and
alternative / post-change one. It must be strictly positive for ST, SR and CU.
Currently, GLRCU does not support the minimum gap, and this param will be ignored.}

\item{delta_upper}{Maximum gap between null / pre-change space and
alternative / post-change one. It must be strictly positive for ST, SR and CU.
Currently, GLRCU does not support the maximum gap, and this param will be ignored.}
}
\value{
A list of
\enumerate{
\item Boolean indicating whether it is acceptable or not.
\item Character describing why it is not acceptable.
\item Updated delta_upper for the case where the original input was NULL
}
}
\description{
For each method and family, check whether delta parameters are within
expected range with respect to the pre-change parameter.
}
