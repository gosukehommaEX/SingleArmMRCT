# Regional Consistency Probability for Single-Arm MRCT (Time-to-Event Endpoint)

Calculate the regional consistency probability (RCP) for time-to-event
endpoints using the hazard ratio (HR) in single-arm multi-regional
clinical trials (MRCTs) using the Effect Retention Approach (ERA).

Event times are modelled by the exponential distribution with hazard
rate \\\lambda\\ (treatment) relative to a known historical control
hazard \\\lambda_0\\. The treatment effect is expressed as a hazard
ratio: \\HR = \lambda / \lambda_0 \< 1\\ (benefit). Two effect scales
are considered:

- Log-HR scale: \\\log(\widehat{HR}\_j) = \log(\hat{\lambda}\_j /
  \lambda_0)\\.

- Linear-HR scale: \\1 - \widehat{HR}\_j = 1 - \hat{\lambda}\_j /
  \lambda_0\\.

Two evaluation methods are supported (for each scale):

- Method 1: Effect retention approach. Evaluates whether Region 1
  retains at least a fraction PI of the overall treatment effect.
  Log-HR: \\\log(\widehat{HR}\_1) \< \pi \times \log(\widehat{HR})\\;
  Linear-HR: \\(1 - \widehat{HR}\_1) \> \pi \times (1 - \widehat{HR})\\.

- Method 2: Simultaneous benefit across all regions. Evaluates whether
  all regional hazard ratios are below 1: \\\widehat{HR}\_j \< 1\\ for
  all \\j\\. (Equivalent for both log-HR and linear-HR scales.)

Two calculation approaches are available:

- `"formula"`: Closed-form solution based on normal approximation for
  \\\log(\widehat{HR})\\ and the delta method for the linear-HR scale
  (Hayashi and Itoh 2018). Method 1 uses a two-block decomposition
  (Region 1 vs regions 2..J combined), which is valid for \\J \geq 2\\.
  Method 2 supports \\J \geq 2\\ regions.

- `"simulation"`: Monte Carlo simulation using individual patient data
  with person-years estimation of the hazard rate. Supports \\J \geq 2\\
  regions.

## Usage

``` r
rcp1armHazardRatio(
  lambda,
  lambda0,
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
  Must be positive. Under exponential distribution, median survival =
  \\\log(2) / \lambda\\.

- lambda0:

  Numeric scalar. Known hazard rate for the historical control (null
  hypothesis reference value). Must be positive.

- Nj:

  Integer vector. Sample sizes for each region. For example, `c(10, 90)`
  indicates Region 1 has 10 subjects and Region 2 has 90 subjects. All
  elements must be positive integers.

- t_a:

  Numeric scalar. Accrual period (patient enrollment duration). Must be
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
  closed-form solution or `"simulation"` for Monte Carlo simulation.
  Default is `"formula"`.

- nsim:

  Positive integer. Number of Monte Carlo iterations. Used only when
  `approach = "simulation"`. Default is `10000`.

- seed:

  Non-negative integer. Random seed for reproducibility. Used only when
  `approach = "simulation"`. Default is `1`.

## Value

An object of class `"rcp1armHazardRatio"`, which is a list containing:

- `approach`:

  Calculation approach used (`"formula"` or `"simulation"`).

- `nsim`:

  Number of Monte Carlo iterations (`NULL` for `"formula"` approach).

- `lambda`:

  True hazard rate under the alternative hypothesis.

- `lambda0`:

  Historical control hazard rate.

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

- `Method1_logHR`:

  RCP using Method 1 (log-HR scale).

- `Method1_linearHR`:

  RCP using Method 1 (linear-HR scale).

- `Method2`:

  RCP using Method 2 (all regions show benefit; identical for log-HR and
  linear-HR scales).

## References

Hayashi R, Itoh Y (2018). A reexamination of Japanese sample size
calculation for multiregional clinical trial evaluating survival
endpoint. *Pharmaceutical Statistics*, 17(1): 46–55.

Wu J (2015). Sample size calculation for the one-sample log-rank test.
*Pharmaceutical Statistics*, 14(1): 26–33.

## Examples

``` r
# Example 1: Closed-form solution with N = 100, Region 1 has 10 subjects
result1 <- rcp1armHazardRatio(
  lambda         = log(2) / 10,
  lambda0        = log(2) / 5,
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
#> Endpoint : Time-to-Event (Hazard Ratio)
#> 
#>    Approach       : Closed-Form Solution
#>    True Hazard    : lambda  = 0.069315
#>    Control Hazard : lambda0 = 0.138629
#>    Sample Size    : Nj      = (10, 90)
#>    Total Size     : N       = 100
#>    Accrual Period : t_a     = 3.00
#>    Follow-up      : t_f     = 10.00
#>    Study Duration : tau     = 13.00
#>    Dropout Hazard : lambda_d = NA
#>    Threshold      : PI      = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall):
#>       Log-HR based    : 0.8007
#>       Linear-HR based : 0.8358
#>    Method 2 (All Regions Show Benefit):
#>       HR < 1          : 0.9478
#> 

# Example 2: Monte Carlo simulation with N = 100, Region 1 has 10 subjects
result2 <- rcp1armHazardRatio(
  lambda         = log(2) / 10,
  lambda0        = log(2) / 5,
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
#> Endpoint : Time-to-Event (Hazard Ratio)
#> 
#>    Approach       : Simulation-Based (nsim = 10000)
#>    True Hazard    : lambda  = 0.069315
#>    Control Hazard : lambda0 = 0.138629
#>    Sample Size    : Nj      = (10, 90)
#>    Total Size     : N       = 100
#>    Accrual Period : t_a     = 3.00
#>    Follow-up      : t_f     = 10.00
#>    Study Duration : tau     = 13.00
#>    Dropout Hazard : lambda_d = NA
#>    Threshold      : PI      = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall):
#>       Log-HR based    : 0.8067
#>       Linear-HR based : 0.8444
#>    Method 2 (All Regions Show Benefit):
#>       HR < 1          : 0.9562
#> 
```
