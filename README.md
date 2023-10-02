
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stcpR6: Sequential Test and Change-Point detection algorithms based on E-values / E-detectors

<!-- badges: start -->

[![R-CMD-check](https://github.com/shinjaehyeok/stcpR6/workflows/R-CMD-check/badge.svg)](https://github.com/shinjaehyeok/stcpR6/actions)
<!-- badges: end -->

stcpR6 is an R package built to run nonparametric sequential tests and
online change point detection algorithms in [SRR
21’](https://arxiv.org/abs/2010.08082) and [SRR
23’](https://arxiv.org/abs/2203.03532). This package supports APIs of
nonparametric sequential test and online change-point detection for
streams of univariate sub-Gaussian, binary, and bounded random
variables.

Disclaimer: This R package is a personal project and is not affiliated
with my company. It is provided “as is” without warranty of any kind,
express or implied. The author does not assume any liability for any
errors or omissions in this package.

## Installation

You can install the development version of stcpR6 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("shinjaehyeok/stcpR6")
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this.  -->
