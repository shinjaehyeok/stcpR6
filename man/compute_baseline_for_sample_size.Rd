% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_baseline.R
\name{compute_baseline_for_sample_size}
\alias{compute_baseline_for_sample_size}
\title{Compute baseline parameters given target variance process bounds.}
\usage{
compute_baseline_for_sample_size(
  alpha,
  v_upper,
  v_lower,
  psi_fn_list = generate_sub_G_fn(),
  skip_g_alpha = TRUE,
  v_min = 1,
  k_max = 200,
  tol = 1e-10
)
}
\arguments{
\item{alpha}{ARL parameter in (0,1)}

\item{v_upper}{Upper bound of the target variance process bound}

\item{v_lower}{Lower bound of the target variance process bound.}

\item{psi_fn_list}{A list of R functions that computes psi and psi_star functions. Can be generated by \code{generate_sub_G_fn()} or counterparts for sub_B and sub_E.}

\item{skip_g_alpha}{If true, we do not compute g_alpha and use log(1/alpha) instead.}

\item{v_min}{A lower bound of v function in the baseline process. Default is \code{1}.}

\item{k_max}{Positive integer to determine the maximum number of baselines. Default is \code{200}.}

\item{tol}{Tolerance of root-finding, positive numeric. Default is 1e-10.}
}
\value{
A list of 1. Parameters of baseline processes, 2. Mixing weights, 3. Auxiliary values for computation.
}
\description{
Given target variance process bounds for confidence sequences, compute baseline parameters.
}
