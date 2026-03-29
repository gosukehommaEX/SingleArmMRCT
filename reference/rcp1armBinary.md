# Regional Consistency Probability for Single-Arm MRCT (Binary Endpoint)

Calculate the regional consistency probability (RCP) for binary
endpoints in single-arm multi-regional clinical trials (MRCTs) using the
Effect Retention Approach (ERA).

Two evaluation methods are supported:

- Method 1: Effect retention approach. Evaluates whether Region 1
  retains at least a fraction PI of the overall treatment effect:
  \\Pr\[(\hat{p}\_1 - p_0) \> \pi \times (\hat{p} - p_0)\]\\.

- Method 2: Simultaneous positivity across all regions. Evaluates
  whether all regional response rates exceed the null value:
  \\Pr\[\hat{p}\_j \> p_0 \text{ for all } j\]\\.

Two calculation approaches are available:

- `"formula"`: Exact closed-form solution via full enumeration of the
  binomial joint distribution. Method 1 uses a two-block decomposition
  (Region 1 vs regions 2..J combined), which is valid for \\J \geq 2\\.
  Method 2 supports \\J \geq 2\\ regions.

- `"simulation"`: Monte Carlo simulation. Supports \\J \geq 2\\ regions.

## Usage

``` r
rcp1armBinary(
  p,
  p0,
  Nj,
  PI = 0.5,
  approach = "formula",
  nsim = 10000,
  seed = 1
)
```

## Arguments

- p:

  Numeric scalar. True response rate under the alternative hypothesis.
  Must be in \\(0, 1)\\.

- p0:

  Numeric scalar. Null hypothesis response rate (baseline or historical
  control). Must be in \\\[0, 1)\\.

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

An object of class `"rcp1armBinary"`, which is a list containing:

- `approach`:

  Calculation approach used (`"formula"` or `"simulation"`).

- `nsim`:

  Number of Monte Carlo iterations (`NULL` for `"formula"` approach).

- `p`:

  True response rate under the alternative hypothesis.

- `p0`:

  Null hypothesis response rate.

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
# Example 1: Exact solution with N = 100, Region 1 has 10 subjects
result1 <- rcp1armBinary(
  p  = 0.5,
  p0 = 0.2,
  Nj = c(10, 90),
  PI = 0.5,
  approach = "formula"
)
print(result1)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Binary
#> 
#>    Approach      : Exact Solution
#>    Response Rate : p  = 0.5000
#>    Null Rate     : p0 = 0.2000
#>    Sample Size   : Nj = (10, 90)
#>    Total Size    : N  = 100
#>    Threshold     : PI = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall) : 0.8309
#>    Method 2 (All Regions > p0)    : 0.9453
#> 

# Example 2: Monte Carlo simulation with N = 100, Region 1 has 10 subjects
result2 <- rcp1armBinary(
  p    = 0.5,
  p0   = 0.2,
  Nj   = c(10, 90),
  PI   = 0.5,
  approach = "simulation",
  nsim = 10000,
  seed = 1
)
print(result2)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Binary
#> 
#>    Approach      : Simulation-Based (nsim = 10000)
#>    Response Rate : p  = 0.5000
#>    Null Rate     : p0 = 0.2000
#>    Sample Size   : Nj = (10, 90)
#>    Total Size    : N  = 100
#>    Threshold     : PI = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall) : 0.8277
#>    Method 2 (All Regions > p0)    : 0.9427
#> 
```
