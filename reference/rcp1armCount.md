# Regional Consistency Probability for Single-Arm MRCT (Count Endpoint)

Calculate the regional consistency probability (RCP) for count
(overdispersed) endpoints in single-arm multi-regional clinical trials
(MRCTs) using the Effect Retention Approach (ERA).

Count data are modelled by the negative binomial distribution, and the
treatment effect is expressed as a rate ratio (RR) relative to a
historical control rate \\\lambda_0\\. Two effect scales are considered:

- Log-RR scale: \\\log(\widehat{RR}\_j) = \log(\hat{\lambda}\_j /
  \lambda_0)\\.

- Linear-RR scale: \\1 - \widehat{RR}\_j = 1 - \hat{\lambda}\_j /
  \lambda_0\\.

Two evaluation methods are supported (for each scale):

- Method 1: Effect retention approach. Evaluates whether Region 1
  retains at least a fraction PI of the overall treatment effect.
  Log-RR: \\\log(\widehat{RR}\_1) \< \pi \times \log(\widehat{RR})\\;
  Linear-RR: \\(1 - \widehat{RR}\_1) \> \pi \times (1 - \widehat{RR})\\.

- Method 2: Simultaneous benefit across all regions. Evaluates whether
  all regional rate ratios are below 1: \\\widehat{RR}\_j \< 1\\ for all
  \\j\\. (Equivalent for both log-RR and linear-RR scales.)

Two calculation approaches are available:

- `"formula"`: Exact closed-form solution via full enumeration of the
  negative binomial joint distribution. Method 1 uses a two-block
  decomposition (Region 1 vs regions 2..J combined), which is valid for
  \\J \geq 2\\. Method 2 supports \\J \geq 2\\ regions.

- `"simulation"`: Monte Carlo simulation. Supports \\J \geq 2\\ regions.

## Usage

``` r
rcp1armCount(
  lambda,
  lambda0,
  dispersion,
  Nj,
  PI = 0.5,
  approach = "formula",
  nsim = 10000,
  seed = 1
)
```

## Arguments

- lambda:

  Numeric scalar. Expected count per patient under the alternative
  hypothesis. Must be positive.

- lambda0:

  Numeric scalar. Expected count per patient under the historical
  control (null hypothesis reference value). Must be positive.

- dispersion:

  Numeric scalar. Dispersion parameter (size) of the negative binomial
  distribution, assumed common across all regions. Smaller values
  indicate greater overdispersion. Must be positive.

- Nj:

  Integer vector. Sample sizes for each region. For example, `c(10, 90)`
  indicates Region 1 has 10 subjects and Region 2 has 90 subjects. All
  elements must be positive integers.

- PI:

  Numeric scalar. Prespecified effect retention threshold for Method 1.
  Typically \\\pi \geq 0.5\\. Must be in \\\[0, 1\]\\. Default is `0.5`.

- approach:

  Character scalar. Calculation approach: `"formula"` for the exact
  solution or `"simulation"` for Monte Carlo simulation. Default is
  `"formula"`.

- nsim:

  Positive integer. Number of Monte Carlo iterations. Used only when
  `approach = "simulation"`. Default is `10000`.

- seed:

  Non-negative integer. Random seed for reproducibility. Used only when
  `approach = "simulation"`. Default is `1`.

## Value

An object of class `"rcp1armCount"`, which is a list containing:

- `approach`:

  Calculation approach used (`"formula"` or `"simulation"`).

- `nsim`:

  Number of Monte Carlo iterations (`NULL` for `"formula"` approach).

- `lambda`:

  Expected count per patient under the alternative hypothesis.

- `lambda0`:

  Expected count per patient under the historical control.

- `dispersion`:

  Dispersion parameter.

- `Nj`:

  Sample sizes for each region.

- `PI`:

  Effect retention threshold.

- `Method1_logRR`:

  RCP using Method 1 (log-RR scale).

- `Method1_linearRR`:

  RCP using Method 1 (linear-RR scale).

- `Method2`:

  RCP using Method 2 (all regions show benefit; identical for log-RR and
  linear-RR scales).

## Examples

``` r
# Example 1: Exact solution with N = 100, Region 1 has 10 subjects
result1 <- rcp1armCount(
  lambda     = 2,
  lambda0    = 3,
  dispersion = 1,
  Nj         = c(10, 90),
  PI         = 0.5,
  approach   = "formula"
)
print(result1)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Count (Negative Binomial)
#> 
#>    Approach       : Exact Solution
#>    Expected Count : lambda     = 2.000000
#>    Control Count  : lambda0    = 3.000000
#>    Dispersion     : dispersion = 1.000000
#>    Sample Size    : Nj         = (10, 90)
#>    Total Size     : N          = 100
#>    Threshold      : PI         = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall):
#>       Log-RR based    : 0.7481
#>       Linear-RR based : 0.7675
#>    Method 2 (All Regions Show Benefit):
#>       RR < 1          : 0.8845
#> 

# Example 2: Monte Carlo simulation with N = 100, Region 1 has 10 subjects
result2 <- rcp1armCount(
  lambda     = 2,
  lambda0    = 3,
  dispersion = 1,
  Nj         = c(10, 90),
  PI         = 0.5,
  approach   = "simulation",
  nsim       = 10000,
  seed       = 1
)
print(result2)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Count (Negative Binomial)
#> 
#>    Approach       : Simulation-Based (nsim = 10000)
#>    Expected Count : lambda     = 2.000000
#>    Control Count  : lambda0    = 3.000000
#>    Dispersion     : dispersion = 1.000000
#>    Sample Size    : Nj         = (10, 90)
#>    Total Size     : N          = 100
#>    Threshold      : PI         = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall):
#>       Log-RR based    : 0.7448
#>       Linear-RR based : 0.7636
#>    Method 2 (All Regions Show Benefit):
#>       RR < 1          : 0.8810
#> 
```
