# SingleArmMRCT 0.1.0

## Initial release

### New functions

**Calculation functions** — each supports `approach = "formula"` (closed-form or semi-analytical) and `approach = "simulation"` (Monte Carlo), and returns both Method 1 and Method 2 RCP:

- `rcp1armContinuous()`: Continuous endpoint. Method 1 uses normal approximation via two-block decomposition; Method 2 uses the product of independent regional probabilities.
- `rcp1armBinary()`: Binary endpoint. Method 1 uses exact enumeration of the binomial joint distribution; Method 2 uses the product of independent binomial tail probabilities.
- `rcp1armCount()`: Count endpoint (negative binomial model). Provides RCP on both the log-RR scale and the linear-RR scale via the reproducibility property of the negative binomial distribution.
- `rcp1armHazardRatio()`: Time-to-event endpoint (hazard ratio). Closed-form solution based on normal approximation for log(HR) and the delta method for the linear-HR scale (Hayashi and Itoh 2018). Supports accrual period, follow-up period, and exponential dropout.
- `rcp1armMilestoneSurvival()`: Milestone survival endpoint. Variance of the Kaplan-Meier estimator derived from Greenwood's formula (Tang 2022). Returns `"closed-form"` when `t_eval <= t_f` and `"numerical-integration"` otherwise.
- `rcp1armRMST()`: Restricted mean survival time (RMST) endpoint. Variance formula via integration of the squared survival difference. Returns `"closed-form"` when `tau_star <= t_f` and `"numerical-integration"` otherwise.

**Plot functions** — each generates a faceted ggplot2 object overlaying formula and simulation results for Method 1 and Method 2 as a function of the regional allocation proportion f₁:

- `plot_rcp1armContinuous()`
- `plot_rcp1armBinary()`
- `plot_rcp1armCount()`
- `plot_rcp1armHazardRatio()`
- `plot_rcp1armMilestoneSurvival()`
- `plot_rcp1armRMST()`

**S3 print methods** — `print()` methods for all six result classes:

- `print.rcp1armContinuous()`
- `print.rcp1armBinary()`
- `print.rcp1armCount()`
- `print.rcp1armHazardRatio()`
- `print.rcp1armMilestoneSurvival()`
- `print.rcp1armRMST()`
