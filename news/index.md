# Changelog

## SingleArmMRCT 0.1.1

CRAN release: 2026-04-01

### Resubmission to CRAN

- Added `Language: en-GB` to DESCRIPTION.
- Expanded abbreviations (RCP, MRCT, RMST, MHLW) in DESCRIPTION.
- Fixed invalid URL in DESCRIPTION (BugReports and URL fields now point
  to the correct GitHub repository: gosukehommaEX/SingleArmMRCT).
- Wrapped long-running examples for
  [`plot_rcp1armRMST()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/plot_rcp1armRMST.md)
  and
  [`plot_rcp1armMilestoneSurvival()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/plot_rcp1armMilestoneSurvival.md)
  in `\donttest{}`.

------------------------------------------------------------------------

## SingleArmMRCT 0.1.0

### Initial release

#### New functions

**Calculation functions** — each supports `approach = "formula"`
(closed-form or semi-analytical) and `approach = "simulation"` (Monte
Carlo), and returns both Method 1 and Method 2 RCP:

- [`rcp1armContinuous()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/rcp1armContinuous.md):
  Continuous endpoint. Method 1 uses normal approximation via two-block
  decomposition; Method 2 uses the product of independent regional
  probabilities.
- [`rcp1armBinary()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/rcp1armBinary.md):
  Binary endpoint. Method 1 uses exact enumeration of the binomial joint
  distribution; Method 2 uses the product of independent binomial tail
  probabilities.
- [`rcp1armCount()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/rcp1armCount.md):
  Count endpoint (negative binomial model). Provides RCP on both the
  log-RR scale and the linear-RR scale via the reproducibility property
  of the negative binomial distribution.
- [`rcp1armHazardRatio()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/rcp1armHazardRatio.md):
  Time-to-event endpoint (hazard ratio). Closed-form solution based on
  normal approximation for log(HR) and the delta method for the
  linear-HR scale (Hayashi and Itoh 2018). Supports accrual period,
  follow-up period, and exponential dropout.
- [`rcp1armMilestoneSurvival()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/rcp1armMilestoneSurvival.md):
  Milestone survival endpoint. Variance of the Kaplan-Meier estimator
  derived from Greenwood’s formula (Tang 2022). Returns `"closed-form"`
  when `t_eval <= t_f` and `"numerical-integration"` otherwise.
- [`rcp1armRMST()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/rcp1armRMST.md):
  Restricted mean survival time (RMST) endpoint. Variance formula via
  integration of the squared survival difference. Returns
  `"closed-form"` when `tau_star <= t_f` and `"numerical-integration"`
  otherwise.

**Plot functions** — each generates a faceted ggplot2 object overlaying
formula and simulation results for Method 1 and Method 2 as a function
of the regional allocation proportion f₁:

- [`plot_rcp1armContinuous()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/plot_rcp1armContinuous.md)
- [`plot_rcp1armBinary()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/plot_rcp1armBinary.md)
- [`plot_rcp1armCount()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/plot_rcp1armCount.md)
- [`plot_rcp1armHazardRatio()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/plot_rcp1armHazardRatio.md)
- [`plot_rcp1armMilestoneSurvival()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/plot_rcp1armMilestoneSurvival.md)
- [`plot_rcp1armRMST()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/plot_rcp1armRMST.md)

**S3 print methods** — [`print()`](https://rdrr.io/r/base/print.html)
methods for all six result classes:

- [`print.rcp1armContinuous()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/print.rcp1armContinuous.md)
- [`print.rcp1armBinary()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/print.rcp1armBinary.md)
- [`print.rcp1armCount()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/print.rcp1armCount.md)
- [`print.rcp1armHazardRatio()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/print.rcp1armHazardRatio.md)
- [`print.rcp1armMilestoneSurvival()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/print.rcp1armMilestoneSurvival.md)
- [`print.rcp1armRMST()`](https://gosukehommaEX.github.io/SingleArmMRCT/reference/print.rcp1armRMST.md)
