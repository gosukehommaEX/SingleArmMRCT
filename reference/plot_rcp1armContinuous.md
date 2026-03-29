# Plot Regional Consistency Probability for Single-Arm MRCT (Continuous Endpoint)

Generate a faceted plot of Regional Consistency Probability (RCP) as a
function of the regional allocation proportion \\f_1\\ for continuous
endpoints. Formula and simulation results are shown together for both
Method 1 and Method 2. Facet columns correspond to total sample sizes
specified in `N_vec`.

Regional sample sizes are allocated as: \\N\_{j1} = \lfloor N \times f_1
\rfloor\\ and \\N\_{j2} = \cdots = N\_{jJ} = (N - N\_{j1}) / (J - 1)\\.

## Usage

``` r
plot_rcp1armContinuous(
  mu = 0.5,
  mu0 = 0.1,
  sd = 1,
  PI = 0.5,
  N_vec = c(20, 40, 100),
  J = 3,
  f1_seq = seq(0.1, 0.9, by = 0.1),
  nsim = 10000,
  seed = 1,
  base_size = 28
)
```

## Arguments

- mu:

  Numeric scalar. True mean under the alternative hypothesis. Default is
  `0.5`.

- mu0:

  Numeric scalar. Null hypothesis mean. Default is `0.1`.

- sd:

  Numeric scalar. True standard deviation. Must be positive. Default is
  `1`.

- PI:

  Numeric scalar. Effect retention threshold for Method 1. Must be in
  \\\[0, 1\]\\. Default is `0.5`.

- N_vec:

  Integer vector. Total sample sizes for each facet column. Default is
  `c(20, 40, 100)`.

- J:

  Positive integer (\>= 2). Number of regions. Default is `3`.

- f1_seq:

  Numeric vector. Sequence of Region 1 allocation proportions. Each
  value must be in \\(0, 1)\\. Default is `seq(0.1, 0.9, by = 0.1)`.

- nsim:

  Positive integer. Number of Monte Carlo iterations for simulation.
  Default is `10000`.

- seed:

  Non-negative integer. Random seed for simulation. Default is `1`.

- base_size:

  Positive numeric. Base font size in points passed to
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html). Use
  larger values (e.g., `28`) for presentation slides and smaller values
  (e.g., `11`) for vignettes or reports. Default is `28`.

## Value

A ggplot2 object.

## Examples

``` r
p <- plot_rcp1armContinuous(
  mu  = 0.5,
  mu0 = 0.1,
  sd  = 1,
  PI  = 0.5,
  N_vec = c(20, 40, 100),
  J   = 3
)
print(p)

```
