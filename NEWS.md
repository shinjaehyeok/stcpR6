# stcpR6 (development version)

# stcpR6 0.9.8
* Two sided e-values are now based on two one-sided e-evalues tuned for alpha/2.
* Sample size based NormalCS tuning is now based on log(1/alpha).
* We can still tune it based on g_alpha by setting skip_g_alpha = FALSE.

# stcpR6 0.9.7

* Improve test coverage for GLR-CUSUM methods
* Refactor Stcp constructor
* Simplify input check logic

# stcpR6 0.9.6

* Start to track updates via NEWS.md
* Reflect minor Cpp header file updates
* Add example for NormalCS class
* Minor patch on NormalCS class implementation
