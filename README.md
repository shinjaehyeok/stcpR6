
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

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this.  -->

## Installation

You can install the development version of stcpR6 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("shinjaehyeok/stcpR6")
```

## Example

### Setup: tracking prediction accuracy.

- Under pre-change, $Y_t = t + ϵ_t$ for $t = 1, 2, \dots, \nu$.
- Under post-change, $Y_t = t + 0.01 * (t-\nu)^2 + ϵ_t$ for
  $t = \nu +1, \nu+2, \dots$.
- Prediction: Simple linear regression based on $(t, Y_t)$ history.

### Build E-detector on residuals.

We will use E-detector for bounded random variables with the following
setup

- Pre-change:
  $\mathbb{E}\left\[Y_t -\hat{Y}_t | Y_1, \dots, Y_{t-1}\right\] \leq 0$
- Post-change :
  $\mathbb{E}\left\[Y_t -\hat{Y}_t | Y_1, \dots, Y_{t-1}\right\] > 0$.
- Minimum practically interesting gap between pre- and post-means :
  $\Delta = 1$.
- Average run length (ARL) level : ARL $\geq 1000$. i.e., False alert
  once in 1000 times.
- Cap residuals into $[-20, 20]$.

``` r
kCap <- 20 
kARL <- 1000

normalize_obs <- function(y) {
  obs <- (y + kCap) / (2 * kCap)
  obs[obs > 1] <- 1
  obs[obs < 0] <- 0
  return(obs)
}

e_detector <- stcpR6::Stcp$new(
    method = "SR",
    family = "Bounded",
    alternative = "greater",
    threshold = log(kARL),
    m_pre = normalize_obs(0),
    delta_lower = 1 / kCap / 2
  )
```

### Case 1. No change happened (Normal condition)

``` r
# Set random seed for reproducibility
set.seed(100)
# No change
sig <- 10
v <- 200

y_history <- seq(1,v) + rnorm(v, 0, sig)

# Simple regression predictor
predict_y <- function(i) {
  if (i <= 2) return (0)
  j <- i-1
  coeff <- lm.fit(x = cbind(rep(1, j),1:j), y = y_history[1:j])$coefficients
  return(coeff[1] + coeff[2] * i)
}

y_hat <- sapply(seq_along(y_history), predict_y)

plot(seq_along(y_history), y_history, type = "l", 
xlab = "t", ylab = "Y_t", main = "Sample history (Black) and Prediction (Red)")
lines(seq_along(y_hat), y_hat, col = 2)
```

<img src="man/figures/README-fig1-1.png" width="100%" />

``` r

plot(seq_along(y_history), y_history-y_hat, type = "l", 
xlab = "t", ylab = "Residual", main = "Prediction error")
```

<img src="man/figures/README-fig1-2.png" width="100%" />

Apply E-detector on normalized residuals

``` r
res <- normalize_obs(as.numeric(y_history-y_hat))
# In real applications, 
# we should use e_detector$updateLogValues(res_i) to update e-deetector
# for each newly observed residual res_i.
# In this example, however, we use a vectorized update for simplicity.

# Initialize the model.
e_detector$reset()
# Compute the log values of e-detectors. 
log_values <- e_detector$updateAndReturnHistories(res) 

# Print a summary of updates
print(e_detector)
#> stcp Model:
#> - Method:  SR 
#> - Family:  Bounded 
#> - Alternative:  greater 
#> - Alpha:  0.001 
#> - m_pre:  0.5 
#> - Num. of mixing components:  169 
#> - Obs. have been passed:  200 
#> - Current log value:  5.068833 
#> - Is stopped before:  FALSE 
#> - Stopped time:  0

# Plot the log values and stopping threshold
plot(seq_along(log_values), log_values, type = "l",
xlab = "t", ylab = "log value", ylim = c(0, 3 * e_detector$getThreshold()))
abline(h = e_detector$getThreshold(), col = 2)
```

<img src="man/figures/README-fig2-1.png" width="100%" />

Note that the log values of e-detector have remained below the
threshold. E-detector didn’t trigger any alert.

### Case 2. Change happened at t = 100

``` r
# Set random seed for reproducibility
set.seed(100)
# Change point: 100
v <- 100

# Pre-change history
y_pre <- seq(1,v) + rnorm(v, 0, sig)

# Post-change history
y_post <- seq(v+1, 2*v) + 0.01 * seq(1,v)^2 + rnorm(v, 0, sig)
y_history <- c(y_pre, y_post)

# Sample history
y_history <- c(y_pre, y_post)

y_hat <- sapply(seq_along(y_history), predict_y)

plot(seq_along(y_history), y_history, type = "l", 
xlab = "t", ylab = "Y_t", main = "Sample history (Black) and Prediction (Red)")
lines(seq_along(y_hat), y_hat, col = 2)
abline(v = v, lty = 2)
```

<img src="man/figures/README-fig3-1.png" width="100%" />

``` r

plot(seq_along(y_history), y_history-y_hat, type = "l", 
xlab = "t", ylab = "Residual", main = "Prediction error")
abline(v = v, lty = 2)
```

<img src="man/figures/README-fig3-2.png" width="100%" />

Apply E-detector on normalized residuals

``` r
res <- normalize_obs(as.numeric(y_history-y_hat))
# In real applications, 
# we should use e_detector$updateLogValues(res_i) to update e-deetector
# for each newly observed residual res_i.
# In this example, however, we use a vectorized update for simplicity.

# Initialize the model.
e_detector$reset()
# Compute the log values of e-detectors. 
log_values <- e_detector$updateAndReturnHistories(res) 

# Print a summary of updates
print(e_detector)
#> stcp Model:
#> - Method:  SR 
#> - Family:  Bounded 
#> - Alternative:  greater 
#> - Alpha:  0.001 
#> - m_pre:  0.5 
#> - Num. of mixing components:  169 
#> - Obs. have been passed:  200 
#> - Current log value:  45.61752 
#> - Is stopped before:  TRUE 
#> - Stopped time:  136

# Plot the log values and stopping threshold
plot(seq_along(log_values), log_values, type = "l",
xlab = "t", ylab = "log value", ylim = c(0, 3 * e_detector$getThreshold()))
abline(h = e_detector$getThreshold() , col = 2)
abline(v = v, lty = 2)
abline(v = e_detector$getStoppedTime(), col = 2, lty = 2)
```

<img src="man/figures/README-fig4-1.png" width="100%" />

In this case, the log values have crossed the threshold at time 136 and
triggered alert. Detection delay was 136 - 100 = 36.

Note that triggering alert at time 136 would be non-trivial if we only
check the residual trend as below.

``` r
plot(seq_along(y_history), y_history-y_hat, type = "l", 
xlab = "t", ylab = "Residual", main = "Prediction error")
abline(v = v, lty = 2)
abline(v = e_detector$getStoppedTime(), col = 2, lty = 2)
```

<img src="man/figures/README-fig5-1.png" width="100%" />
