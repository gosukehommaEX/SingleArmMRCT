# Regional Consistency Probability for Single-Arm MRCT (RMST Endpoint)

Calculate the regional consistency probability (RCP) for restricted mean
survival time (RMST) endpoints in single-arm multi-regional clinical
trials (MRCTs) using the Effect Retention Approach (ERA).

Event times are modelled by the exponential distribution with hazard
rate \\\lambda\\ (treatment). The treatment effect is expressed as the
difference in RMST at a prespecified truncation time \\\tau^\*\\:
\\\delta = \mu(\tau^\*) - \mu_0(\tau^\*)\\, where \\\mu_0(\tau^\*)\\ is
the historical control RMST at \\\tau^\*\\.

The regional RMST estimator is defined as the area under the
Kaplan-Meier curve up to \\\tau^\*\\: \\\hat{\mu}\_j(\tau^\*) =
\int_0^{\tau^\*} \hat{S}\_j(t)\\dt\\.

Two evaluation methods are supported:

- Method 1: Effect retention approach. Evaluates whether Region 1
  retains at least a fraction PI of the overall treatment effect:
  \\Pr\[(\hat{\mu}\_1 - \mu_0) \> \pi \times (\hat{\mu} - \mu_0)\]\\.

- Method 2: Simultaneous benefit across all regions. Evaluates whether
  all regional RMST estimates exceed the historical control RMST:
  \\Pr\[\hat{\mu}\_j \> \mu_0 \text{ for all } j\]\\.

Two calculation approaches are available:

- `"formula"`: Normal approximation based on exponential event times and
  exponential dropout. When \\\tau^\* \leq t_f\\ (the most common
  setting in practice), the administrative censoring survival function
  equals 1 on \\\[0, \tau^\*\]\\ and the variance integral has a
  closed-form solution. When \\\tau^\* \> t_f\\, the integral is
  evaluated numerically via
  [`integrate`](https://rdrr.io/r/stats/integrate.html). Method 1 uses a
  two-block decomposition (Region 1 vs regions 2..J combined), valid for
  \\J \geq 2\\. Method 2 supports \\J \geq 2\\ regions.

- `"simulation"`: Monte Carlo simulation. The RMST for each simulation
  replicate is computed as the exact area under the Kaplan-Meier step
  function (sum of rectangles). Supports \\J \geq 2\\ regions.

## Usage

``` r
rcp1armRMST(
  lambda,
  tau_star,
  mu0,
  Nj,
  t_a,
  t_f,
  lambda_dropout = NULL,
  PI = 0.5,
  approach = "formula",
  nsim = 10000,
  seed = 1
)
```

## Arguments

- lambda:

  Numeric scalar. True hazard rate under the alternative hypothesis.
  Must be positive.

- tau_star:

  Numeric scalar. Truncation time for RMST calculation. Must be positive
  and no greater than the total study duration \\\tau = t_a + t_f\\.

- mu0:

  Numeric scalar. Historical control RMST at `tau_star`, i.e.,
  \\\mu_0(\tau^\*) = \int_0^{\tau^\*} S_0(t)\\dt\\. Must be positive and
  less than `tau_star`.

- Nj:

  Integer vector. Sample sizes for each region. For example, `c(10, 90)`
  indicates Region 1 has 10 subjects and Region 2 has 90 subjects. All
  elements must be positive integers.

- t_a:

  Numeric scalar. Accrual period (patient enrolment duration). Must be
  positive.

- t_f:

  Numeric scalar. Follow-up period (additional follow-up after accrual
  ends). Must be positive.

- lambda_dropout:

  Numeric scalar or `NULL`. Dropout hazard rate. If `NULL` (default), no
  dropout is assumed. If specified, dropout times follow an exponential
  distribution with rate `lambda_dropout`.

- PI:

  Numeric scalar. Prespecified effect retention threshold for Method 1.
  Typically \\\pi \geq 0.5\\. Must be in \\\[0, 1\]\\. Default is `0.5`.

- approach:

  Character scalar. Calculation approach: `"formula"` for the normal
  approximation or `"simulation"` for Monte Carlo simulation. Default is
  `"formula"`.

- nsim:

  Positive integer. Number of Monte Carlo iterations. Used only when
  `approach = "simulation"`. Default is `10000`.

- seed:

  Non-negative integer. Random seed for reproducibility. Used only when
  `approach = "simulation"`. Default is `1`.

## Value

An object of class `"rcp1armRMST"`, which is a list containing:

- `approach`:

  Calculation approach used (`"formula"` or `"simulation"`).

- `formula_type`:

  For `approach = "formula"`: either `"closed-form"` (when \\\tau^\*
  \leq t_f\\) or `"numerical-integration"` (when \\\tau^\* \> t_f\\).
  `NULL` for simulation.

- `nsim`:

  Number of Monte Carlo iterations (`NULL` for `"formula"` approach).

- `lambda`:

  True hazard rate under the alternative hypothesis.

- `Nj`:

  Sample sizes for each region.

- `t_a`:

  Accrual period.

- `t_f`:

  Follow-up period.

- `tau`:

  Total study duration (\\\tau = t_a + t_f\\).

- `lambda_dropout`:

  Dropout hazard rate (`NA` if `NULL`).

- `PI`:

  Effect retention threshold.

- `tau_star`:

  Truncation time for RMST.

- `mu0`:

  Historical control RMST at `tau_star`.

- `mu_est`:

  True RMST under the alternative: \\\mu(\tau^\*) = (1 - e^{-\lambda
  \tau^\*}) / \lambda\\.

- `Method1`:

  RCP using Method 1 (effect retention).

- `Method2`:

  RCP using Method 2 (all regions positive).

## Details

**Variance formula.** The asymptotic variance of the RMST estimator
\\\hat{\mu}\_j(\tau^\*)\\ for a single arm with \\n\\ subjects is: \$\$
\mathrm{Var}(\hat{\mu}\_j(\tau^\*)) \approx \frac{1}{n} \int_0^{\tau^\*}
\frac{\left(e^{-\lambda t} - e^{-\lambda \tau^\*}\right)^2} {\lambda
\cdot e^{-\lambda t} \cdot G(t)} \\dt, \$\$ where \\G(t) =
P(\mathrm{observed\\ time} \geq t)\\ is the overall censoring survival
function: \$\$ G(t) = e^{-\lambda_d t} \times G_a(t), \quad G_a(t) =
\begin{cases} 1 & 0 \leq t \leq t_f, \\ (\tau-t)/t_a & t_f \< t \leq
\tau. \end{cases} \$\$ **Closed-form variance (\\\tau^\* \leq t_f\\).**
When \\\tau^\* \leq t_f\\, \\G_a(t) = 1\\ on \\\[0, \tau^\*\]\\ and the
integrand reduces to \\e^{\lambda_d t}(1 - e^{-\lambda(\tau^\*-t)})^2 /
\lambda\\. Expanding and integrating term by term: \$\$ v(\tau^\*) =
\frac{1}{\lambda}\Bigl\[ A(\lambda_d) - 2 e^{-\lambda\tau^\*}
A(\lambda_d + \lambda) + e^{-2\lambda\tau^\*} A(\lambda_d + 2\lambda)
\Bigr\], \$\$ where \\A(r) = (e^{r\tau^\*} - 1)/r\\ for \\r \> 0\\ and
\\A(0) = \tau^\*\\. When there is no dropout (\\\lambda_d = 0\\) this
further reduces to: \$\$ v(\tau^\*) = \tau^\* -
\frac{2(1-e^{-\lambda\tau^\*})}{\lambda} +
\frac{1-e^{-2\lambda\tau^\*}}{2\lambda}. \$\$

## References

Wu J (2015). Sample size calculation for the one-sample log-rank test.
*Pharmaceutical Statistics*, 14(1): 26–33.

## Examples

``` r
# Example 1: Closed-form solution (tau_star <= t_f) with N = 100
lam   <- log(2) / 10
lam0  <- log(2) / 5
tstar <- 8
mu0_val <- (1 - exp(-lam0 * tstar)) / lam0

result1 <- rcp1armRMST(
  lambda         = lam,
  tau_star       = tstar,
  mu0            = mu0_val,
  Nj             = c(10, 90),
  t_a            = 3,
  t_f            = 10,
  lambda_dropout = NULL,
  PI             = 0.5,
  approach       = "formula"
)
print(result1)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Restricted Mean Survival Time (RMST)
#> 
#>    Approach       : Closed-Form Solution
#>    True Hazard    : lambda  = 0.069315
#>    Sample Size    : Nj      = (10, 90)
#>    Total Size     : N       = 100
#>    Accrual Period : t_a     = 3.00
#>    Follow-up      : t_f     = 10.00
#>    Study Duration : tau     = 13.00
#>    Dropout Hazard : lambda_d = NA
#>    Threshold      : PI      = 0.5000
#>    Trunc. Time    : tau*    = 8.00
#>    Control RMST   : mu0     = 4.8339
#>    True RMST      : mu_est  = 6.1408
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall)  : 0.7767
#>    Method 2 (All Regions > mu0)    : 0.9284
#> 

# Example 2: Monte Carlo simulation with N = 100
result2 <- rcp1armRMST(
  lambda         = lam,
  tau_star       = tstar,
  mu0            = mu0_val,
  Nj             = c(10, 90),
  t_a            = 3,
  t_f            = 10,
  lambda_dropout = NULL,
  PI             = 0.5,
  approach       = "simulation",
  nsim           = 10000,
  seed           = 1
)
print(result2)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Restricted Mean Survival Time (RMST)
#> 
#>    Approach       : Simulation-Based (nsim = 10000)
#>    True Hazard    : lambda  = 0.069315
#>    Sample Size    : Nj      = (10, 90)
#>    Total Size     : N       = 100
#>    Accrual Period : t_a     = 3.00
#>    Follow-up      : t_f     = 10.00
#>    Study Duration : tau     = 13.00
#>    Dropout Hazard : lambda_d = NA
#>    Threshold      : PI      = 0.5000
#>    Trunc. Time    : tau*    = 8.00
#>    Control RMST   : mu0     = 4.8339
#>    True RMST      : mu_est  = 6.1408
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall)  : 0.7920
#>    Method 2 (All Regions > mu0)    : 0.9335
#> 
```
