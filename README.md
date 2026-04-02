# SingleArmMRCT <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
![R](https://img.shields.io/badge/R-%3E%3D4.1.0-blue)
[![CRAN status](https://www.r-pkg.org/badges/version/SingleArmMRCT)](https://CRAN.R-project.org/package=SingleArmMRCT)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/SingleArmMRCT)](https://cranlogs.r-pkg.org/badges/grand-total/SingleArmMRCT)
[![CRAN downloads monthly](https://cranlogs.r-pkg.org/badges/SingleArmMRCT)](https://cranlogs.r-pkg.org/badges/SingleArmMRCT)
![License](https://img.shields.io/badge/license-MIT-yellow)
<!-- badges: end -->

## Overview

**SingleArmMRCT** provides functions to calculate and visualise the **Regional Consistency Probability (RCP)** for single-arm multi-regional clinical trials (MRCTs) using the **Effect Retention Approach (ERA)**.

The package addresses a critical methodological gap: current Japanese MHLW Method 1 and Method 2 consistency criteria were originally developed for two-arm trials, yet single-arm trials increasingly form the basis of regulatory submissions, particularly in oncology. This package extends classical approaches to the single-arm setting across six endpoint types.

### Supported endpoints

| Endpoint type | Calculation function | Plot function |
|---|---|---|
| Continuous | `rcp1armContinuous()` | `plot_rcp1armContinuous()` |
| Binary | `rcp1armBinary()` | `plot_rcp1armBinary()` |
| Count (negative binomial) | `rcp1armCount()` | `plot_rcp1armCount()` |
| Time-to-event (hazard ratio) | `rcp1armHazardRatio()` | `plot_rcp1armHazardRatio()` |
| Milestone survival | `rcp1armMilestoneSurvival()` | `plot_rcp1armMilestoneSurvival()` |
| Restricted mean survival time (RMST) | `rcp1armRMST()` | `plot_rcp1armRMST()` |

### Consistency evaluation methods

- **Method 1** (Effect Retention): Evaluates whether Region 1 retains at least a fraction pi of the overall treatment effect.
- **Method 2** (Simultaneous Positivity): Evaluates whether all regional estimates exceed the null value simultaneously.

### Calculation approaches

Each function supports two approaches:

- **`"formula"`**: Closed-form or semi-analytical solution based on normal approximation.
- **`"simulation"`**: Monte Carlo simulation.

---

## Installation

```r
install.packages("SingleArmMRCT")
```

---

## Quick start

### Continuous endpoint

```r
library(SingleArmMRCT)

# Closed-form solution: N = 100, Region 1 has 10 subjects (f1 = 0.1)
result <- rcp1armContinuous(
  mu  = 0.5,
  mu0 = 0.1,
  sd  = 1,
  Nj  = c(10, 90),
  PI  = 0.5,
  approach = "formula"
)
print(result)
```

### Binary endpoint

```r
result <- rcp1armBinary(
  p  = 0.5,
  p0 = 0.2,
  Nj = c(10, 90),
  PI = 0.5,
  approach = "formula"
)
print(result)
```

### Time-to-event endpoint (hazard ratio)

```r
result <- rcp1armHazardRatio(
  lambda         = log(2) / 10,
  lambda0        = log(2) / 5,
  Nj             = c(10, 90),
  t_a            = 3,
  t_f            = 10,
  lambda_dropout = NULL,
  PI             = 0.5,
  approach       = "formula"
)
print(result)
```

### RMST endpoint

```r
lam0    <- log(2) / 5
tstar   <- 8
mu0_val <- (1 - exp(-lam0 * tstar)) / lam0

result <- rcp1armRMST(
  lambda   = log(2) / 10,
  tau_star = tstar,
  mu0      = mu0_val,
  Nj       = c(10, 90),
  t_a      = 3,
  t_f      = 10,
  PI       = 0.5,
  approach = "formula"
)
print(result)
```

---

## Visualisation

Each endpoint has a corresponding plot function that generates a faceted plot of RCP as a function of the regional allocation proportion f₁, overlaying formula and simulation results for both Method 1 and Method 2.

```r
# Continuous endpoint: RCP vs f1 for N = 20, 40, 100 with J = 3 regions
p <- plot_rcp1armContinuous(
  mu    = 0.5,
  mu0   = 0.1,
  sd    = 1,
  PI    = 0.5,
  N_vec = c(20, 40, 100),
  J     = 3
)
print(p)
```

```r
# Milestone survival endpoint
p <- plot_rcp1armMilestoneSurvival(
  lambda = log(2) / 10,
  t_eval = 8,
  S0     = exp(-log(2) * 8 / 5),
  t_a    = 3,
  t_f    = 10,
  PI     = 0.5,
  N_vec  = c(20, 40, 100),
  J      = 3
)
print(p)
```

---

## Parameter conventions

| Symbol | Meaning | Notes |
|---|---|---|
| `Nj` | Integer vector of regional sample sizes | e.g., `c(10, 90)` for J = 2 regions |
| `PI` | Effect retention threshold pi | Typically ≥ 0.5; default 0.5 |
| `f1` | Regional allocation proportion of Region 1 | f₁ = Nj[1] / sum(Nj) |
| `t_a` | Accrual period | Time-to-event endpoints only |
| `t_f` | Follow-up period | Time-to-event endpoints only |
| `lambda_dropout` | Dropout hazard rate | `NULL` = no dropout |

---

## References

Hayashi N, Itoh Y (2017). A re-examination of Japanese sample size calculation for multi-regional clinical trial evaluating survival endpoint. *Japanese Journal of Biometrics*, 38(2): 79--92. <https://doi.org/10.5691/jjb.38.79>

Homma G (2024). Cautionary note on regional consistency evaluation in multiregional clinical trials with binary outcomes. *Pharmaceutical Statistics*, 23(3):385--398. <https://doi.org/10.1002/pst.2358>

Wu J (2015). Sample size calculation for the one-sample log-rank test. *Pharmaceutical Statistics*, 14(1): 26--33. <https://doi.org/10.1002/pst.1654>

---

## License

MIT © 2025 Gosuke Homma
