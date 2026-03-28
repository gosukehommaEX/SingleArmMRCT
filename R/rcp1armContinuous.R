#' Regional Consistency Probability for Single-Arm MRCT (Continuous Endpoint)
#'
#' @description
#' Calculate the regional consistency probability (RCP) for continuous endpoints
#' in single-arm multi-regional clinical trials (MRCTs) using the Effect Retention
#' Approach (ERA).
#'
#' Two evaluation methods are supported:
#' \itemize{
#'   \item Method 1: Effect retention approach. Evaluates whether Region 1 retains
#'     at least a fraction PI of the overall treatment effect:
#'     \eqn{Pr[(\hat{\mu}_1 - \mu_0) > \pi \times (\hat{\mu} - \mu_0)]}.
#'   \item Method 2: Simultaneous positivity across all regions.
#'     Evaluates whether all regional estimates exceed the null value:
#'     \eqn{Pr[\hat{\mu}_j > \mu_0 \text{ for all } j]}.
#' }
#'
#' Two calculation approaches are available:
#' \itemize{
#'   \item \code{"formula"}: Closed-form analytical solution based on normal
#'     approximation. Method 1 uses a two-block decomposition (Region 1 vs
#'     regions 2..J combined), which is valid for \eqn{J \geq 2}.
#'     Method 2 supports \eqn{J \geq 2} regions.
#'   \item \code{"simulation"}: Monte Carlo simulation. Supports \eqn{J \geq 2} regions.
#' }
#'
#' @param mu Numeric scalar. True mean under the alternative hypothesis.
#' @param mu0 Numeric scalar. Null hypothesis mean (baseline or historical control).
#'   The treatment effect is defined as \eqn{\delta = \mu - \mu_0}.
#' @param sd Numeric scalar. True standard deviation, assumed common across all
#'   regions. Must be positive.
#' @param Nj Integer vector. Sample sizes for each region. For example,
#'   \code{c(10, 90)} indicates Region 1 has 10 subjects and Region 2 has 90
#'   subjects. All elements must be positive integers.
#' @param PI Numeric scalar. Prespecified effect retention threshold for Method 1.
#'   Typically \eqn{\pi \geq 0.5}. Must be in \eqn{[0, 1]}. Default is \code{0.5}.
#' @param approach Character scalar. Calculation approach: \code{"formula"} for the
#'   closed-form solution or \code{"simulation"} for Monte Carlo simulation.
#'   Default is \code{"formula"}.
#' @param nsim Positive integer. Number of Monte Carlo iterations. Used only when
#'   \code{approach = "simulation"}. Default is \code{10000}.
#' @param seed Non-negative integer. Random seed for reproducibility. Used only
#'   when \code{approach = "simulation"}. Default is \code{1}.
#'
#' @return An object of class \code{"rcp1armContinuous"}, which is a list containing:
#' \describe{
#'   \item{\code{approach}}{Calculation approach used (\code{"formula"} or
#'     \code{"simulation"}).}
#'   \item{\code{nsim}}{Number of Monte Carlo iterations (\code{NULL} for
#'     \code{"formula"} approach).}
#'   \item{\code{mu}}{True mean under the alternative hypothesis.}
#'   \item{\code{mu0}}{Null hypothesis mean.}
#'   \item{\code{sd}}{Standard deviation.}
#'   \item{\code{Nj}}{Sample sizes for each region.}
#'   \item{\code{PI}}{Effect retention threshold.}
#'   \item{\code{Method1}}{RCP using Method 1 (effect retention).}
#'   \item{\code{Method2}}{RCP using Method 2 (all regions positive).}
#' }
#'
#' @importFrom stats rnorm pnorm
#'
#' @examples
#' # Example 1: Closed-form solution with N = 100, Region 1 has 10 subjects
#' result1 <- rcp1armContinuous(
#'   mu  = 0.5,
#'   mu0 = 0.1,
#'   sd  = 1,
#'   Nj  = c(10, 90),
#'   PI  = 0.5,
#'   approach = "formula"
#' )
#' print(result1)
#'
#' # Example 2: Monte Carlo simulation with N = 100, Region 1 has 10 subjects
#' result2 <- rcp1armContinuous(
#'   mu   = 0.5,
#'   mu0  = 0.1,
#'   sd   = 1,
#'   Nj   = c(10, 90),
#'   PI   = 0.5,
#'   approach = "simulation",
#'   nsim = 10000,
#'   seed = 1
#' )
#' print(result2)
#'
#' @export
rcp1armContinuous <- function(mu,
                              mu0,
                              sd,
                              Nj,
                              PI       = 0.5,
                              approach = "formula",
                              nsim     = 1e4,
                              seed     = 1) {

  # ========== Input Validation ==========
  if (!is.numeric(mu) || length(mu) != 1) {
    stop("mu must be a single numeric value")
  }
  if (!is.numeric(mu0) || length(mu0) != 1) {
    stop("mu0 must be a single numeric value")
  }
  if (!is.numeric(sd) || length(sd) != 1 || sd <= 0) {
    stop("sd must be a single positive number")
  }
  if (!is.numeric(Nj) || any(Nj <= 0) || any(Nj != as.integer(Nj))) {
    stop("Nj must be a vector of positive integers")
  }
  if (!is.numeric(PI) || length(PI) != 1 || PI < 0 || PI > 1) {
    stop("PI must be a single numeric value in [0, 1]")
  }
  if (!approach %in% c("formula", "simulation")) {
    stop('approach must be either "formula" or "simulation"')
  }
  if (approach == "simulation") {
    if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 || nsim %% 1 != 0) {
      stop("nsim must be a single positive integer")
    }
    if (!is.numeric(seed) || length(seed) != 1 || seed < 0 || seed %% 1 != 0) {
      stop("seed must be a single non-negative integer")
    }
  }

  # ========== Common Setup ==========
  N     <- sum(Nj)
  J     <- length(Nj)
  fj    <- Nj / N
  delta <- mu - mu0

  # ========== Calculation ==========
  if (approach == "formula") {
    # Method 1: Pr[(hat_mu1 - mu0) > PI * (hat_mu - mu0)] = Pr(D > 0)
    # D = (1 - PI*f1) * hat_delta1 - PI*(1-f1) * hat_delta_others
    # where hat_delta_others is the pooled estimate for regions 2..J combined.
    # This two-block decomposition (Region 1 vs all others) is valid for J >= 2.
    # Under homogeneity: E[D] = (1 - PI) * delta
    # Var[D] = (1 - PI*f1)^2 * sigma^2/N1 + (PI*(1-f1))^2 * sigma^2/(N-N1)
    mean_d <- (1 - PI) * delta
    var_d  <- (1 - PI * fj[1])^2 * (sd^2 / Nj[1]) +
      (PI * (1 - fj[1]))^2 * (sd^2 / (N - Nj[1]))
    sd_d   <- sqrt(var_d)

    # RCP: Pr(D > 0) = Phi(E[D] / SD[D])
    Method1 <- stats::pnorm(mean_d / sd_d)

    # Method 2: Pr(hat_muj > mu0 for all j) = prod_j Pr(hat_muj > mu0)
    # Regions are independent, so the product holds for J >= 2.
    # hat_muj ~ N(mu, sd^2 / Nj), so Pr(hat_muj > mu0) = Phi(delta * sqrt(Nj) / sd)
    Method2 <- prod(stats::pnorm(delta * sqrt(Nj) / sd))

    nsim_out <- NULL

  } else {
    # ----- Monte Carlo Simulation -----
    set.seed(seed)

    # Generate regional sample means: hat_muj ~ N(mu, sd^2 / Nj[j])
    # Dimensions: nsim rows x J columns
    hat_muj_mat <- matrix(NA_real_, nrow = nsim, ncol = J)
    for (j in seq_len(J)) {
      hat_muj_mat[, j] <- stats::rnorm(nsim, mean = mu, sd = sd / sqrt(Nj[j]))
    }

    # Overall sample mean: hat_mu = sum(Nj * hat_muj) / N
    hat_mu <- as.vector(hat_muj_mat %*% Nj) / N

    # Method 1: Pr[(hat_mu1 - mu0) > PI * (hat_mu - mu0)]
    Method1 <- mean((hat_muj_mat[, 1] - mu0) > PI * (hat_mu - mu0))

    # Method 2: Pr[hat_muj > mu0 for all j]
    Method2 <- mean(rowSums(hat_muj_mat > mu0) == J)

    nsim_out <- nsim
  }

  # ========== Output ==========
  result <- list(
    approach = approach,
    nsim     = nsim_out,
    mu       = mu,
    mu0      = mu0,
    sd       = sd,
    Nj       = Nj,
    PI       = PI,
    Method1  = Method1,
    Method2  = Method2
  )
  class(result) <- "rcp1armContinuous"
  return(result)
}


#' Print Method for rcp1armContinuous Objects
#'
#' @param x An object of class \code{"rcp1armContinuous"}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
print.rcp1armContinuous <- function(x, ...) {
  cat("\nRegional Consistency Probability for Single-Arm MRCT\n")
  cat("Endpoint : Continuous\n\n")

  if (x$approach == "simulation") {
    cat(sprintf("   Approach    : Simulation-Based (nsim = %d)\n", x$nsim))
  } else {
    cat("   Approach    : Closed-Form Solution\n")
  }

  cat(sprintf("   Target Mean : mu  = %.4f\n", x$mu))
  cat(sprintf("   Null Mean   : mu0 = %.4f\n", x$mu0))
  cat(sprintf("   Std. Dev.   : sd  = %.4f\n", x$sd))
  cat(sprintf("   Sample Size : Nj  = (%s)\n", paste(x$Nj, collapse = ", ")))
  cat(sprintf("   Total Size  : N   = %d\n",   sum(x$Nj)))
  cat(sprintf("   Threshold   : PI  = %.4f\n", x$PI))

  cat("\nConsistency Probabilities:\n")
  cat(sprintf("   Method 1 (Region 1 vs Overall)  : %.4f\n", x$Method1))
  cat(sprintf("   Method 2 (All Regions > mu0)    : %.4f\n", x$Method2))
  cat("\n")

  invisible(x)
}
