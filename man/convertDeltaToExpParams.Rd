% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{convertDeltaToExpParams}
\alias{convertDeltaToExpParams}
\title{converted input deltas to parameters for exponential baselines}
\usage{
convertDeltaToExpParams(
  family,
  alternative,
  threshold,
  m_pre,
  delta_lower,
  delta_upper,
  k_max
)
}
\arguments{
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

\item{threshold}{Stopping threshold. We recommend to use log(1/alpha)
for "ST" and "SR" methods where alpha is a testing level or 1/ARL.
for "CU" and "GRLCU", we recommend to tune the threshold by using
domain-specific sampler to hit the target ARL.}

\item{m_pre}{The boundary of mean parameter in null / pre-change space}

\item{delta_lower}{Minimum gap between null / pre-change space and
alternative / post-change one. It must be strictly positive for ST, SR and CU.
Currently, GLRCU does not support the minimum gap, and this param will be ignored.}

\item{delta_upper}{Maximum gap between null / pre-change space and
alternative / post-change one. It must be strictly positive for ST, SR and CU.
Currently, GLRCU does not support the maximum gap, and this param will be ignored.}

\item{k_max}{Positive integer to determine the maximum number of baselines.
For GLRCU method, it is used as the lookup window size for GLRCU statistics.}
}
\value{
A list of weights and lambdas
}
\description{
For each exponential baseline family, convert delta range into
corresponding lambdas and weights.
}
