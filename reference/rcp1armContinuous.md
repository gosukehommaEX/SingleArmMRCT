# Regional Consistency Probability for Single-Arm MRCT (Continuous Endpoint)

Calculate the regional consistency probability (RCP) for continuous
endpoints in single-arm multi-regional clinical trials (MRCTs) using the
Effect Retention Approach (ERA).

Two evaluation methods are supported:

- Method 1: Effect retention approach. Evaluates whether Region 1
  retains at least a fraction PI of the overall treatment effect:
  \\Pr\[(\hat{\mu}\_1 - \mu_0) \> \pi \times (\hat{\mu} - \mu_0)\]\\.

- Method 2: Simultaneous positivity across all regions. Evaluates
  whether all regional estimates exceed the null value:
  \\Pr\[\hat{\mu}\_j \> \mu_0 \text{ for all } j\]\\.

Two calculation approaches are available:

- `"formula"`: Closed-form analytical solution based on normal
  approximation. Method 1 uses a two-block decomposition (Region 1 vs
  regions 2..J combined), which is valid for \\J \geq 2\\. Method 2
  supports \\J \geq 2\\ regions.

- `"simulation"`: Monte Carlo simulation. Supports \\J \geq 2\\ regions.

## Usage

``` r
rcp1armContinuous(
  mu,
  mu0,
  sd,
  Nj,
  PI = 0.5,
  approach = "formula",
  nsim = 10000,
  seed = 1
)
```

## Arguments

- mu:

  Numeric scalar. True mean under the alternative hypothesis.

- mu0:

  Numeric scalar. Null hypothesis mean (baseline or historical control).
  The treatment effect is defined as \\\delta = \mu - \mu_0\\.

- sd:

  Numeric scalar. True standard deviation, assumed common across all
  regions. Must be positive.

- Nj:

  Integer vector. Sample sizes for each region. For example, `c(10, 90)`
  indicates Region 1 has 10 subjects and Region 2 has 90 subjects. All
  elements must be positive integers.

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

An object of class `"rcp1armContinuous"`, which is a list containing:

- `approach`:

  Calculation approach used (`"formula"` or `"simulation"`).

- `nsim`:

  Number of Monte Carlo iterations (`NULL` for `"formula"` approach).

- `mu`:

  True mean under the alternative hypothesis.

- `mu0`:

  Null hypothesis mean.

- `sd`:

  Standard deviation.

- `Nj`:

  Sample sizes for each region.

- `PI`:

  Effect retention threshold.

- `Method1`:

  RCP using Method 1 (effect retention).

- `Method2`:

  RCP using Method 2 (all regions positive).

## Examples

``` r
# Example 1: Closed-form solution with N = 100, Region 1 has 10 subjects
result1 <- rcp1armContinuous(
  mu  = 0.5,
  mu0 = 0.1,
  sd  = 1,
  Nj  = c(10, 90),
  PI  = 0.5,
  approach = "formula"
)
print(result1)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Continuous
#> 
#>    Approach    : Closed-Form Solution
#>    Target Mean : mu  = 0.5000
#>    Null Mean   : mu0 = 0.1000
#>    Std. Dev.   : sd  = 1.0000
#>    Sample Size : Nj  = (10, 90)
#>    Total Size  : N   = 100
#>    Threshold   : PI  = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall)  : 0.7446
#>    Method 2 (All Regions > mu0)    : 0.8970
#> 

# Example 2: Monte Carlo simulation with N = 100, Region 1 has 10 subjects
result2 <- rcp1armContinuous(
  mu   = 0.5,
  mu0  = 0.1,
  sd   = 1,
  Nj   = c(10, 90),
  PI   = 0.5,
  approach = "simulation",
  nsim = 10000,
  seed = 1
)
print(result2)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Continuous
#> 
#>    Approach    : Simulation-Based (nsim = 10000)
#>    Target Mean : mu  = 0.5000
#>    Null Mean   : mu0 = 0.1000
#>    Std. Dev.   : sd  = 1.0000
#>    Sample Size : Nj  = (10, 90)
#>    Total Size  : N   = 100
#>    Threshold   : PI  = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall)  : 0.7421
#>    Method 2 (All Regions > mu0)    : 0.8922
#> 
```
