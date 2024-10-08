## Resubmission
This is a resubmission. In this version I have:

* Two sided e-values are now based on two one-sided e-evalues tuned for alpha/2.
* Sample size based NormalCS tuning is now based on log(1/alpha).
* We can still tune it based on g_alpha by setting skip_g_alpha = FALSE.

* Simplify input check logic

## R CMD check results

0 errors | 0 warnings | 0 note

