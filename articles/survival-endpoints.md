# Survival Endpoints: Hazard Ratio, Milestone Survival, and RMST

This vignette describes Regional Consistency Probability (RCP)
calculations for three survival endpoint types: **hazard ratio**,
**milestone survival probability**, and **restricted mean survival time
(RMST)**. All three endpoints share a common trial design framework, and
event times are modelled by the exponential distribution.

------------------------------------------------------------------------

## Common trial design framework

For all survival endpoints, the following parameters define the trial
design.

| Parameter        | Symbol                 | Description                                                                        |
|------------------|------------------------|------------------------------------------------------------------------------------|
| `t_a`            | $t_{a}$                | Accrual period; patients enrol uniformly over $\left\lbrack 0,t_{a} \right\rbrack$ |
| `t_f`            | $t_{f}$                | Follow-up period after accrual closes                                              |
| `tau`            | $\tau = t_{a} + t_{f}$ | Total study duration (computed internally)                                         |
| `lambda`         | $\lambda$              | True hazard rate under the alternative (exponential model)                         |
| `lambda0`        | $\lambda_{0}$          | Historical control hazard rate                                                     |
| `lambda_dropout` | $\lambda_{d}$          | Dropout hazard rate; `NULL` assumes no dropout (default)                           |

The common parameters below are used throughout the examples in this
vignette.

``` r
lambda  <- log(2) / 10   # treatment arm: median survival = 10
lambda0 <- log(2) / 5    # historical control: median survival = 5
t_a     <- 3             # accrual period
t_f     <- 10            # follow-up period
# True HR = lambda / lambda0 = 0.5
```

------------------------------------------------------------------------

## 1. Hazard Ratio Endpoint

### Statistical model

Under the exponential model with uniform accrual over
$\left\lbrack 0,t_{a} \right\rbrack$ and administrative censoring at
$\tau$, the expected event probability per patient is:

$$\phi = \frac{\lambda}{\lambda + \lambda_{d}}\left\lbrack 1 - \frac{e^{- {(\lambda + \lambda_{d})}t_{f}} - e^{- {(\lambda + \lambda_{d})}\tau}}{\left( \lambda + \lambda_{d} \right)\, t_{a}} \right\rbrack$$

The expected number of events in Region $j$ is $E_{j} = N_{j}\,\phi$,
and the log-hazard ratio estimator has the approximate distribution:

$$\log\left( {\widehat{HR}}_{j} \right) \sim N\!\left( \log(HR),\;\frac{1}{E_{j}} \right)$$

### Consistency criteria

**Method 1 (log-HR scale):**

Letting $\delta = \log(HR) < 0$, $E_{1} = N_{1}\phi$, and
$E_{- 1} = \left( N - N_{1} \right)\phi$:

$$\text{RCP}_{1,\log} = \Phi\!\left( \frac{- (1 - \pi)\,\delta}{\sqrt{\left( 1 - \pi f_{1} \right)^{2}/E_{1} + \{\pi\left( 1 - f_{1} \right)\}^{2}/E_{- 1}}} \right)$$

**Method 1 (linear-HR scale)** (Hayashi and Itoh 2018):

$$\text{RCP}_{1,\text{lin}} = \Pr\!\left\lbrack \,\left( 1 - {\widehat{HR}}_{1} \right) \geq \pi\,\left( 1 - \widehat{HR} \right)\, \right\rbrack$$

This is derived via the delta method. Define
$g = \log\left( {\widehat{HR}}_{1} \right) - \log\left( 1 - \pi + \pi\,\widehat{HR} \right)$.
Under homogeneity:

$$E\lbrack g\rbrack = \log(HR) - \log(1 - \pi + \pi\, HR)$$

$${Var}(g) = \frac{1}{E_{\text{total}}}\left\lbrack \frac{\{ 1 - \pi + \pi\left( 1 - f_{1} \right)HR\}^{2}}{f_{1}\,(1 - \pi + \pi\, HR)^{2}} + \frac{\{\pi\left( 1 - f_{1} \right)HR\}^{2}}{\left( 1 - f_{1} \right)\,(1 - \pi + \pi\, HR)^{2}} \right\rbrack$$

and
$\text{RCP}_{1,\text{lin}} = \Phi\left( - E\lbrack g\rbrack/\sqrt{{Var}(g)} \right)$.

**Method 2:**

$$\text{RCP}_{2} = \prod\limits_{j = 1}^{J}\Pr\!\left( {\widehat{HR}}_{j} < 1 \right) = \prod\limits_{j = 1}^{J}\Phi\!\left( - \delta\,\sqrt{E_{j}} \right)$$

### Example

``` r
result_f <- rcp1armHazardRatio(
  lambda         = lambda,
  lambda0        = lambda0,
  Nj             = c(20, 80),
  t_a            = t_a,
  t_f            = t_f,
  lambda_dropout = NULL,
  PI             = 0.5,
  approach       = "formula"
)
print(result_f)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Time-to-Event (Hazard Ratio)
#> 
#>    Approach       : Closed-Form Solution
#>    True Hazard    : lambda  = 0.069315
#>    Control Hazard : lambda0 = 0.138629
#>    Sample Size    : Nj      = (20, 80)
#>    Total Size     : N       = 100
#>    Accrual Period : t_a     = 3.00
#>    Follow-up      : t_f     = 10.00
#>    Study Duration : tau     = 13.00
#>    Dropout Hazard : lambda_d = NA
#>    Threshold      : PI      = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall):
#>       Log-HR based    : 0.8935
#>       Linear-HR based : 0.9228
#>    Method 2 (All Regions Show Benefit):
#>       HR < 1          : 0.9892
```

``` r
result_s <- rcp1armHazardRatio(
  lambda         = lambda,
  lambda0        = lambda0,
  Nj             = c(20, 80),
  t_a            = t_a,
  t_f            = t_f,
  lambda_dropout = NULL,
  PI             = 0.5,
  approach       = "simulation",
  nsim           = 10000,
  seed           = 1
)
print(result_s)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Time-to-Event (Hazard Ratio)
#> 
#>    Approach       : Simulation-Based (nsim = 10000)
#>    True Hazard    : lambda  = 0.069315
#>    Control Hazard : lambda0 = 0.138629
#>    Sample Size    : Nj      = (20, 80)
#>    Total Size     : N       = 100
#>    Accrual Period : t_a     = 3.00
#>    Follow-up      : t_f     = 10.00
#>    Study Duration : tau     = 13.00
#>    Dropout Hazard : lambda_d = NA
#>    Threshold      : PI      = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall):
#>       Log-HR based    : 0.9019
#>       Linear-HR based : 0.9320
#>    Method 2 (All Regions Show Benefit):
#>       HR < 1          : 0.9924
```

### Effect of dropout

Specifying `lambda_dropout` reduces the expected event probability
$\phi$, which in turn reduces all RCP values.

``` r
result_dropout <- rcp1armHazardRatio(
  lambda         = lambda,
  lambda0        = lambda0,
  Nj             = c(20, 80),
  t_a            = t_a,
  t_f            = t_f,
  lambda_dropout = 0.05,
  PI             = 0.5,
  approach       = "formula"
)
print(result_dropout)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Time-to-Event (Hazard Ratio)
#> 
#>    Approach       : Closed-Form Solution
#>    True Hazard    : lambda  = 0.069315
#>    Control Hazard : lambda0 = 0.138629
#>    Sample Size    : Nj      = (20, 80)
#>    Total Size     : N       = 100
#>    Accrual Period : t_a     = 3.00
#>    Follow-up      : t_f     = 10.00
#>    Study Duration : tau     = 13.00
#>    Dropout Hazard : lambda_d = 0.050000
#>    Threshold      : PI      = 0.5000
#> 
#> Consistency Probabilities:
#>    Method 1 (Region 1 vs Overall):
#>       Log-HR based    : 0.8656
#>       Linear-HR based : 0.8971
#>    Method 2 (All Regions Show Benefit):
#>       HR < 1          : 0.9793
```

### Visualisation

``` r
plot_rcp1armHazardRatio(
  lambda    = lambda,
  lambda0   = lambda0,
  t_a       = t_a,
  t_f       = t_f,
  PI        = 0.5,
  N_vec     = c(20, 40, 100),
  J         = 3,
  nsim      = 5000,
  seed      = 1,
  base_size = 8
)
```

![Grid plot of RCP versus f1 for a hazard ratio endpoint with HR = 0.5,
showing Method 1 on log-HR and linear-HR scales and Method 2 across N =
20, 40, 100](survival-endpoints_files/figure-html/unnamed-chunk-5-1.png)

------------------------------------------------------------------------

## 2. Milestone Survival Endpoint

### Statistical model

The treatment effect at evaluation time $t_{\text{eval}}$ is:

$$\delta = S\left( t_{\text{eval}} \right) - S_{0} = e^{- \lambda\, t_{\text{eval}}} - S_{0}$$

where $S_{0}$ is the historical control survival rate at
$t_{\text{eval}}$, supplied by the user.

The asymptotic variance of the Kaplan-Meier estimator
${\widehat{S}}_{j}\left( t_{\text{eval}} \right)$ is derived from
Greenwood’s formula:

$${Var}\!\left\lbrack {\widehat{S}}_{j}(t) \right\rbrack \approx \frac{e^{- 2\lambda t}}{N_{j}}\int_{0}^{t}\frac{\lambda\, e^{{(\lambda + \lambda_{d})}u}}{G_{a}(u)}\, du\; = :\;\frac{v_{KM}(t)}{N_{j}}$$

where:

$$G_{a}(u) = \begin{cases}
1 & {0 \leq u \leq t_{f}} \\
{(\tau - u)/t_{a}} & {t_{f} < u \leq \tau}
\end{cases}$$

**Closed-form solution when $t_{\text{eval}} \leq t_{f}$**
($G_{a}(u) = 1$ throughout):

$$v_{KM}(t) = e^{- 2\lambda t} \cdot \frac{\lambda}{\lambda + \lambda_{d}}\left( e^{{(\lambda + \lambda_{d})}\, t} - 1 \right)$$

For $\lambda_{d} = 0$ this reduces to the binomial variance
$S(t)\left( 1 - S(t) \right)$. When $t_{\text{eval}} > t_{f}$, the
integral is evaluated numerically via
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html).

### Consistency criteria

**Method 1:**

$$\text{RCP}_{1} = \Phi\!\left( \frac{(1 - \pi)\,\delta}{\sqrt{\left( 1 - \pi f_{1} \right)^{2}\, v_{KM}/N_{1} + \{\pi\left( 1 - f_{1} \right)\}^{2}\, v_{KM}/\left( N - N_{1} \right)}} \right)$$

**Method 2:**

$$\text{RCP}_{2} = \prod\limits_{j = 1}^{J}\Phi\!\left( \frac{\delta}{\sqrt{v_{KM}/N_{j}}} \right)$$

### Example

``` r
t_eval <- 8
S0     <- exp(-log(2) * t_eval / 5)
cat(sprintf("True S(%g) = %.4f,  S0 = %.4f,  delta = %.4f\n",
            t_eval, exp(-lambda * t_eval), S0,
            exp(-lambda * t_eval) - S0))
#> True S(8) = 0.5743,  S0 = 0.3299,  delta = 0.2445
```

``` r
result_f <- rcp1armMilestoneSurvival(
  lambda         = lambda,
  t_eval         = t_eval,
  S0             = S0,
  Nj             = c(20, 80),
  t_a            = t_a,
  t_f            = t_f,
  lambda_dropout = NULL,
  PI             = 0.5,
  approach       = "formula"
)
print(result_f)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Milestone Survival
#> 
#>    Approach       : Closed-Form Solution (Greenwood)
#>    True Hazard    : lambda  = 0.069315
#>    Sample Size    : Nj      = (20, 80)
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
#>    Method 1 (Region 1 vs Overall)  : 0.8848
#>    Method 2 (All Regions > S0)     : 0.9865
```

Since $t_{\text{eval}} = 8 \leq t_{f} = 10$, the closed-form solution is
used (`formula_type = "closed-form"`). For $t_{\text{eval}} > t_{f}$,
numerical integration is applied automatically.

``` r
result_s <- rcp1armMilestoneSurvival(
  lambda         = lambda,
  t_eval         = t_eval,
  S0             = S0,
  Nj             = c(20, 80),
  t_a            = t_a,
  t_f            = t_f,
  lambda_dropout = NULL,
  PI             = 0.5,
  approach       = "simulation",
  nsim           = 10000,
  seed           = 1
)
print(result_s)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Milestone Survival
#> 
#>    Approach       : Simulation-Based (nsim = 10000)
#>    True Hazard    : lambda  = 0.069315
#>    Sample Size    : Nj      = (20, 80)
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
#>    Method 1 (Region 1 vs Overall)  : 0.8873
#>    Method 2 (All Regions > S0)     : 0.9887
```

### Visualisation

``` r
plot_rcp1armMilestoneSurvival(
  lambda    = lambda,
  t_eval    = t_eval,
  S0        = S0,
  t_a       = t_a,
  t_f       = t_f,
  PI        = 0.5,
  N_vec     = c(20, 40, 100),
  J         = 3,
  nsim      = 5000,
  seed      = 1,
  base_size = 8
)
```

![Line plot of RCP versus f1 for a milestone survival endpoint at t_eval
= 8, showing Method 1 and Method 2 across N = 20, 40,
100](survival-endpoints_files/figure-html/unnamed-chunk-9-1.png)

------------------------------------------------------------------------

## 3. RMST Endpoint

### Statistical model

The restricted mean survival time (RMST) up to truncation time
$\tau^{*}$ is:

$$\mu\left( \tau^{*} \right) = \int_{0}^{\tau^{*}}S(t)\, dt = \frac{1 - e^{- \lambda\tau^{*}}}{\lambda}$$

The asymptotic variance of the regional RMST estimator is:

$${Var}\!\left\lbrack {\widehat{\mu}}_{j}\left( \tau^{*} \right) \right\rbrack \approx \frac{v_{\text{RMST}}\left( \tau^{*} \right)}{N_{j}}$$

where:

$$v_{\text{RMST}}\left( \tau^{*} \right) = \int_{0}^{\tau^{*}}\frac{e^{\lambda_{d}t}\left( 1 - e^{- \lambda{(\tau^{*} - t)}} \right)^{2}}{\lambda\, G_{a}(t)}\, dt$$

**Closed-form solution when $\tau^{*} \leq t_{f}$**, using
$A(r) = \left( e^{r\tau^{*}} - 1 \right)/r$ (and $A(0) = \tau^{*}$):

\$\$ v\_{\text{RMST}}(\tau^\*) = \frac{1}{\lambda}\Bigl\[ A(\lambda_d) -
2\\e^{-\lambda\tau^\*}\\A(\lambda_d + \lambda) +
e^{-2\lambda\tau^\*}\\A(\lambda_d + 2\lambda) \Bigr\] \$\$

For $\lambda_{d} = 0$ this simplifies to:

$$v_{\text{RMST}}\left( \tau^{*} \right) = \tau^{*} - \frac{2\left( 1 - e^{- \lambda\tau^{*}} \right)}{\lambda} + \frac{1 - e^{- 2\lambda\tau^{*}}}{2\lambda}$$

When $\tau^{*} > t_{f}$, the integral is split at $t_{f}$ and evaluated
via [`stats::integrate()`](https://rdrr.io/r/stats/integrate.html).

### Consistency criteria

**Method 1:**

$$\text{RCP}_{1} = \Phi\!\left( \frac{(1 - \pi)\,\left( \mu_{\text{est}} - \mu_{0} \right)}{\sqrt{\left( 1 - \pi f_{1} \right)^{2}\, v/N_{1} + \{\pi\left( 1 - f_{1} \right)\}^{2}\, v/\left( N - N_{1} \right)}} \right)$$

**Method 2:**

$$\text{RCP}_{2} = \prod\limits_{j = 1}^{J}\Phi\!\left( \frac{\mu_{\text{est}} - \mu_{0}}{\sqrt{v/N_{j}}} \right)$$

### Example

``` r
tau_star <- 8
mu0      <- (1 - exp(-lambda0 * tau_star)) / lambda0
mu_est   <- (1 - exp(-lambda  * tau_star)) / lambda
cat(sprintf("True RMST = %.4f,  mu0 = %.4f,  delta = %.4f\n",
            mu_est, mu0, mu_est - mu0))
#> True RMST = 6.1408,  mu0 = 4.8339,  delta = 1.3069
```

``` r
result_f <- rcp1armRMST(
  lambda         = lambda,
  tau_star       = tau_star,
  mu0            = mu0,
  Nj             = c(20, 80),
  t_a            = t_a,
  t_f            = t_f,
  lambda_dropout = NULL,
  PI             = 0.5,
  approach       = "formula"
)
print(result_f)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Restricted Mean Survival Time (RMST)
#> 
#>    Approach       : Closed-Form Solution
#>    True Hazard    : lambda  = 0.069315
#>    Sample Size    : Nj      = (20, 80)
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
#>    Method 1 (Region 1 vs Overall)  : 0.8693
#>    Method 2 (All Regions > mu0)    : 0.9808
```

``` r
result_s <- rcp1armRMST(
  lambda         = lambda,
  tau_star       = tau_star,
  mu0            = mu0,
  Nj             = c(20, 80),
  t_a            = t_a,
  t_f            = t_f,
  lambda_dropout = NULL,
  PI             = 0.5,
  approach       = "simulation",
  nsim           = 10000,
  seed           = 1
)
print(result_s)
#> 
#> Regional Consistency Probability for Single-Arm MRCT
#> Endpoint : Restricted Mean Survival Time (RMST)
#> 
#>    Approach       : Simulation-Based (nsim = 10000)
#>    True Hazard    : lambda  = 0.069315
#>    Sample Size    : Nj      = (20, 80)
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
#>    Method 1 (Region 1 vs Overall)  : 0.8790
#>    Method 2 (All Regions > mu0)    : 0.9809
```

### Visualisation

``` r
plot_rcp1armRMST(
  lambda    = lambda,
  tau_star  = tau_star,
  mu0       = mu0,
  t_a       = t_a,
  t_f       = t_f,
  PI        = 0.5,
  N_vec     = c(20, 40, 100),
  J         = 3,
  nsim      = 5000,
  seed      = 1,
  base_size = 8
)
```

![Line plot of RCP versus f1 for an RMST endpoint with tau_star = 8,
showing Method 1 and Method 2 across N = 20, 40,
100](survival-endpoints_files/figure-html/unnamed-chunk-13-1.png)

------------------------------------------------------------------------

## Summary

| Endpoint           | Effect parameter                                                                                                                             | Benefit direction               | Variance basis                       | Closed-form condition        |
|:-------------------|:---------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------|:-------------------------------------|:-----------------------------|
| Hazard Ratio       | $\log(HR) = \log\left( \lambda/\lambda_{0} \right)$ (Method 1, log-HR scale); $1 - HR = 1 - \lambda/\lambda_{0}$ (Method 1, linear-HR scale) | ${\widehat{HR}}_{j} < 1$        | Expected events via $\phi$ (Wu 2015) | Always                       |
| Milestone Survival | $\delta = e^{- \lambda t_{\text{eval}}} - S_{0}$                                                                                             | ${\widehat{S}}_{j}(t) > S_{0}$  | Greenwood’s formula                  | $t_{\text{eval}} \leq t_{f}$ |
| RMST               | $\delta = \mu\left( \tau^{*} \right) - \mu_{0}\left( \tau^{*} \right)$                                                                       | ${\widehat{\mu}}_{j} > \mu_{0}$ | Squared survival difference integral | $\tau^{*} \leq t_{f}$        |

------------------------------------------------------------------------

## References

Hayashi N, Itoh Y (2017). A re-examination of Japanese sample size
calculation for multi-regional clinical trial evaluating survival
endpoint. *Japanese Journal of Biometrics*, 38(2): 79–92.
<https://doi.org/10.5691/jjb.38.79>

Wu J (2015). Sample size calculation for the one-sample log-rank test.
*Pharmaceutical Statistics*, 14(1): 26–33.
<https://doi.org/10.1002/pst.1654>
