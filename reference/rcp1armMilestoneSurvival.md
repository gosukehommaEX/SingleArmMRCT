# Regional Consistency Probability for Single-Arm MRCT (Milestone Survival Endpoint)

Calculate the regional consistency probability (RCP) for milestone
survival endpoints in single-arm multi-regional clinical trials (MRCTs)
using the Effect Retention Approach (ERA).

Event times are modelled by the exponential distribution with hazard
rate \\\lambda\\ (treatment). The treatment effect at a prespecified
evaluation time \\t\_{\mathrm{eval}}\\ is expressed as the difference in
survival rates: \\\delta = S(t\_{\mathrm{eval}}) -
S_0(t\_{\mathrm{eval}})\\, where \\S_0\\ is the historical control
survival rate at \\t\_{\mathrm{eval}}\\.

Two evaluation methods are supported:

- Method 1: Effect retention approach. Evaluates whether Region 1
  retains at least a fraction PI of the overall treatment effect:
  \\Pr\[(\hat{S}\_1(t) - S_0(t)) \> \pi \times (\hat{S}(t) -
  S_0(t))\]\\.

- Method 2: Simultaneous benefit across all regions. Evaluates whether
  all regional Kaplan-Meier estimates exceed the historical control:
  \\Pr\[\hat{S}\_j(t) \> S_0(t) \text{ for all } j\]\\.

Two calculation approaches are available:

- `"formula"`: Closed-form or semi-analytical solution based on the
  asymptotic variance of the Kaplan-Meier estimator derived from
  Greenwood's formula. When \\t\_{\mathrm{eval}} \leq t_f\\, the
  administrative censoring survival function \\G_a(t) = 1\\ and the
  variance integral has a closed-form solution. When
  \\t\_{\mathrm{eval}} \> t_f\\, the integral is evaluated numerically
  via [`integrate`](https://rdrr.io/r/stats/integrate.html). Method 1
  uses a two-block decomposition (Region 1 vs regions 2..J combined).
  Method 2 supports \\J \geq 2\\ regions.

- `"simulation"`: Monte Carlo simulation with vectorized Kaplan-Meier
  estimation. Supports \\J \geq 2\\ regions.

## Usage

``` r
rcp1armMilestoneSurvival(
  lambda,
  t_eval,
  S0,
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

- t_eval:

  Numeric scalar. Evaluation time point for the milestone survival
  probability. Must be positive.

- S0:

  Numeric scalar. Historical control survival rate at `t_eval`. Must be
  in \\(0, 1\]\\.

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

  Character scalar. Calculation approach: `"formula"` for the
  closed-form (or semi-analytical) solution or `"simulation"` for Monte
  Carlo simulation. Default is `"formula"`.

- nsim:

  Positive integer. Number of Monte Carlo iterations. Used only when
  `approach = "simulation"`. Default is `10000`.

- seed:

  Non-negative integer. Random seed for reproducibility. Used only when
  `approach = "simulation"`. Default is `1`.

## Value

An object of class `"rcp1armMilestoneSurvival"`, which is a list
containing:

- `approach`:

  Calculation approach used (`"formula"` or `"simulation"`).

- `formula_type`:

  For `approach = "formula"`: either `"closed-form"` (when
  \\t\_{\mathrm{eval}} \leq t_f\\) or `"numerical-integration"` (when
  \\t\_{\mathrm{eval}} \> t_f\\). `NULL` for simulation.

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

- `eval_time`:

  Milestone evaluation time point.

- `S0`:

  Historical control survival rate at `eval_time`.

- `S_est`:

  True survival rate under the alternative at `eval_time`: \\S(t) =
  e^{-\lambda t}\\.

- `Method1`:

  RCP using Method 1 (effect retention).

- `Method2`:

  RCP using Method 2 (all regions positive).

## Details

**Variance formula.** The asymptotic variance of the Kaplan-Meier
estimator \\\hat{S}\_j(t)\\ for a single arm with \\N_j\\ subjects is
derived from Greenwood's formula: \$\$ \mathrm{Var}\[\hat{S}\_j(t)\]
\approx \frac{S^2(t)}{N_j} \int_0^t \frac{\lambda(u)}{S(u) \cdot G(u)}
\\ du, \$\$ where \\G(u) = e^{-\lambda_d u} \cdot G_a(u)\\ is the
overall censoring survival function, and \$\$ G_a(u) = \begin{cases} 1 &
0 \leq u \leq t_f, \\ (\tau - u)/t_a & t_f \< u \leq \tau. \end{cases}
\$\$ Under the exponential model \\S(u) = e^{-\lambda u}\\, this
simplifies to: \$\$ \mathrm{Var}\[\hat{S}\_j(t)\] \approx
\frac{e^{-2\lambda t}}{N_j} \int_0^t \frac{\lambda \\ e^{(\lambda +
\lambda_d) u}}{G_a(u)} \\ du. \$\$

**Closed-form solution (\\t\_{\mathrm{eval}} \leq t_f\\).** When
\\t\_{\mathrm{eval}} \leq t_f\\, \\G_a(u) = 1\\ throughout \\\[0,
t\_{\mathrm{eval}}\]\\, and the integral reduces to: \$\$ \int_0^t
\lambda \\ e^{(\lambda + \lambda_d) u} \\ du = \frac{\lambda}{\lambda +
\lambda_d} \bigl( e^{(\lambda + \lambda_d) t} - 1 \bigr). \$\$
Therefore: \$\$ \mathrm{Var}\[\hat{S}\_j(t)\] = \frac{e^{-2\lambda
t}}{N_j} \cdot \frac{\lambda}{\lambda + \lambda_d} \bigl( e^{(\lambda +
\lambda_d) t} - 1 \bigr). \$\$ When there is no dropout (\\\lambda_d =
0\\), this further simplifies to the binomial variance: \$\$
\mathrm{Var}\[\hat{S}\_j(t)\] = \frac{S(t)(1 - S(t))}{N_j}. \$\$

## References

Wu J (2015). Sample size calculation for the one-sample log-rank test.
*Pharmaceutical Statistics*, 14(1): 26–33.

## Examples

``` r
# Example 1: Closed-form solution (t_eval <= t_f) with N = 100
result1 <- rcp1armMilestoneSurvival(
  lambda         = log(2) / 10,
  t_eval         = 8,
  S0             = exp(-log(2) * 8 / 5),
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
#> Endpoint : Milestone Survival
#> 
#>    Approach       : Closed-Form Solution (Greenwood)
#>    True Hazard    : lambda  = 0.069315
#>    Sample Size    : Nj      = (10, 90)
#>    Total Size     : N       = 100
#>    Accrual Period : t_a     = 3.00
#>    Follow-up      : t_f     = 10.00
#>    Study Duration : tau     = 13.00
#>    Dropout Hazard : lambda_d = NA
#>    Threshold      : PI      = 0.5000
#>    Eval Time      : t_eval  = 8.00
#>    Control Surv   : S0      = 0.3299
#>    True Surv      : S_est   = 0.5743
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall)  : 0.7918
#>    Method 2 (All Regions > S0)     : 0.9410
#> 

# Example 2: Monte Carlo simulation with N = 100
result2 <- rcp1armMilestoneSurvival(
  lambda         = log(2) / 10,
  t_eval         = 8,
  S0             = exp(-log(2) * 8 / 5),
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
#> Endpoint : Milestone Survival
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
#>    Eval Time      : t_eval  = 8.00
#>    Control Surv   : S0      = 0.3299
#>    True Surv      : S_est   = 0.5743
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall)  : 0.7922
#>    Method 2 (All Regions > S0)     : 0.9224
#> 
```
