#' Regional Consistency Probability for Single-Arm MRCT (Count Endpoint)
#'
#' @description
#' Calculate the regional consistency probability (RCP) for count (overdispersed)
#' endpoints in single-arm multi-regional clinical trials (MRCTs) using the Effect
#' Retention Approach (ERA).
#'
#' Count data are modelled by the negative binomial distribution, and the treatment
#' effect is expressed as a rate ratio (RR) relative to a historical control rate
#' \eqn{\lambda_0}. Two effect scales are considered:
#' \itemize{
#'   \item Log-RR scale: \eqn{\log(\widehat{RR}_j) = \log(\hat{\lambda}_j / \lambda_0)}.
#'   \item Linear-RR scale: \eqn{1 - \widehat{RR}_j = 1 - \hat{\lambda}_j / \lambda_0}.
#' }
#'
#' Two evaluation methods are supported (for each scale):
#' \itemize{
#'   \item Method 1: Effect retention approach. Evaluates whether Region 1 retains
#'     at least a fraction PI of the overall treatment effect.
#'     Log-RR: \eqn{\log(\widehat{RR}_1) < \pi \times \log(\widehat{RR})};
#'     Linear-RR: \eqn{(1 - \widehat{RR}_1) > \pi \times (1 - \widehat{RR})}.
#'   \item Method 2: Simultaneous benefit across all regions.
#'     Evaluates whether all regional rate ratios are below 1:
#'     \eqn{\widehat{RR}_j < 1} for all \eqn{j}.
#'     (Equivalent for both log-RR and linear-RR scales.)
#' }
#'
#' Two calculation approaches are available:
#' \itemize{
#'   \item \code{"formula"}: Exact closed-form solution via full enumeration of
#'     the negative binomial joint distribution. Method 1 uses a two-block
#'     decomposition (Region 1 vs regions 2..J combined), which is valid for
#'     \eqn{J \geq 2}. Method 2 supports \eqn{J \geq 2} regions.
#'   \item \code{"simulation"}: Monte Carlo simulation. Supports \eqn{J \geq 2} regions.
#' }
#'
#' @param lambda Numeric scalar. Expected count per patient under the alternative
#'   hypothesis. Must be positive.
#' @param lambda0 Numeric scalar. Expected count per patient under the historical
#'   control (null hypothesis reference value). Must be positive.
#' @param dispersion Numeric scalar. Dispersion parameter (size) of the negative
#'   binomial distribution, assumed common across all regions. Smaller values
#'   indicate greater overdispersion. Must be positive.
#' @param Nj Integer vector. Sample sizes for each region. For example,
#'   \code{c(10, 90)} indicates Region 1 has 10 subjects and Region 2 has 90
#'   subjects. All elements must be positive integers.
#' @param PI Numeric scalar. Prespecified effect retention threshold for Method 1.
#'   Typically \eqn{\pi \geq 0.5}. Must be in \eqn{[0, 1]}. Default is \code{0.5}.
#' @param approach Character scalar. Calculation approach: \code{"formula"} for the
#'   exact solution or \code{"simulation"} for Monte Carlo simulation.
#'   Default is \code{"formula"}.
#' @param nsim Positive integer. Number of Monte Carlo iterations. Used only when
#'   \code{approach = "simulation"}. Default is \code{10000}.
#' @param seed Non-negative integer. Random seed for reproducibility. Used only
#'   when \code{approach = "simulation"}. Default is \code{1}.
#'
#' @return An object of class \code{"rcp1armCount"}, which is a list containing:
#' \describe{
#'   \item{\code{approach}}{Calculation approach used (\code{"formula"} or
#'     \code{"simulation"}).}
#'   \item{\code{nsim}}{Number of Monte Carlo iterations (\code{NULL} for
#'     \code{"formula"} approach).}
#'   \item{\code{lambda}}{Expected count per patient under the alternative hypothesis.}
#'   \item{\code{lambda0}}{Expected count per patient under the historical control.}
#'   \item{\code{dispersion}}{Dispersion parameter.}
#'   \item{\code{Nj}}{Sample sizes for each region.}
#'   \item{\code{PI}}{Effect retention threshold.}
#'   \item{\code{Method1_logRR}}{RCP using Method 1 (log-RR scale).}
#'   \item{\code{Method1_linearRR}}{RCP using Method 1 (linear-RR scale).}
#'   \item{\code{Method2}}{RCP using Method 2 (all regions show benefit;
#'     identical for log-RR and linear-RR scales).}
#' }
#'
#' @importFrom stats rnbinom dnbinom qnbinom pnbinom
#'
#' @examples
#' # Example 1: Exact solution with N = 100, Region 1 has 10 subjects
#' result1 <- rcp1armCount(
#'   lambda     = 2,
#'   lambda0    = 3,
#'   dispersion = 1,
#'   Nj         = c(10, 90),
#'   PI         = 0.5,
#'   approach   = "formula"
#' )
#' print(result1)
#'
#' # Example 2: Monte Carlo simulation with N = 100, Region 1 has 10 subjects
#' result2 <- rcp1armCount(
#'   lambda     = 2,
#'   lambda0    = 3,
#'   dispersion = 1,
#'   Nj         = c(10, 90),
#'   PI         = 0.5,
#'   approach   = "simulation",
#'   nsim       = 10000,
#'   seed       = 1
#' )
#' print(result2)
#'
#' @export
rcp1armCount <- function(lambda,
                         lambda0,
                         dispersion,
                         Nj,
                         PI       = 0.5,
                         approach = "formula",
                         nsim     = 1e4,
                         seed     = 1) {

  # ========== Input Validation ==========
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) {
    stop("lambda must be a single positive number")
  }
  if (!is.numeric(lambda0) || length(lambda0) != 1 || lambda0 <= 0) {
    stop("lambda0 must be a single positive number")
  }
  if (!is.numeric(dispersion) || length(dispersion) != 1 || dispersion <= 0) {
    stop("dispersion must be a single positive number")
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
  N  <- sum(Nj)
  J  <- length(Nj)
  fj <- Nj / N

  # ========== Calculation ==========
  if (approach == "formula") {
    # Method 1: Exact enumeration over all (y1, y_others) combinations.
    # y_others represents the pooled count for regions 2..J combined.
    # By the reproducibility of the negative binomial distribution,
    # y_others ~ NB(mu = (N-N1)*lambda, size = (N-N1)*dispersion).
    # This two-block decomposition (Region 1 vs all others) is valid for J >= 2.
    # Method 2: Product of independent regional probabilities, valid for J >= 2.
    #
    # Negative binomial reproducibility property:
    # Sum of Nj i.i.d. NB(mu = lambda, size = dispersion) follows
    # NB(mu = Nj * lambda, size = Nj * dispersion).

    # ----- Region 1 marginal distribution -----
    mu_1   <- Nj[1] * lambda
    size_1 <- Nj[1] * dispersion
    y1_max <- stats::qnbinom(0.9999, mu = mu_1, size = size_1)
    y1_vals <- 0:y1_max
    p_y1    <- stats::dnbinom(y1_vals, mu = mu_1, size = size_1)

    # ----- Regions 2..J combined marginal distribution -----
    n_others    <- N - Nj[1]
    mu_others   <- n_others * lambda
    size_others <- n_others * dispersion
    yo_max  <- stats::qnbinom(0.9999, mu = mu_others, size = size_others)
    yo_vals <- 0:yo_max
    p_yo    <- stats::dnbinom(yo_vals, mu = mu_others, size = size_others)

    # Joint probability matrix and count matrices
    prob_matrix <- outer(p_y1, p_yo)
    y1_matrix   <- matrix(y1_vals, nrow = length(y1_vals), ncol = length(yo_vals))
    yo_matrix   <- matrix(yo_vals, nrow = length(y1_vals), ncol = length(yo_vals),
                          byrow = TRUE)

    # Regional and overall rate estimates per patient
    hat_lambda1_matrix <- y1_matrix / Nj[1]
    hat_lambda_matrix  <- (y1_matrix + yo_matrix) / N

    # Rate ratios relative to historical control
    RR1_matrix    <- hat_lambda1_matrix / lambda0
    RR_all_matrix <- hat_lambda_matrix  / lambda0

    # ----- Method 1: Log-RR -----
    # Condition: log(RR_1) < PI * log(RR_all)
    # Note: RR = 0 (zero counts) yields log(0) = -Inf, which correctly satisfies
    # the condition only when RR_all = 0 as well. To avoid -Inf < PI * (-Inf)
    # being evaluated as FALSE in R (NaN), we handle the degenerate case
    # (y1 = 0 AND y_others = 0) explicitly by treating the condition as FALSE
    # (no observed data to assess consistency).
    both_zero      <- (y1_matrix == 0) & (yo_matrix == 0)
    log_RR1        <- ifelse(RR1_matrix    > 0, log(RR1_matrix),    -Inf)
    log_RR_all     <- ifelse(RR_all_matrix > 0, log(RR_all_matrix), -Inf)
    condition_logRR <- (!both_zero) & (log_RR1 < PI * log_RR_all)
    Method1_logRR  <- sum(prob_matrix * condition_logRR)

    # ----- Method 1: Linear-RR -----
    # Condition: (1 - RR_1) > PI * (1 - RR_all)
    condition_linearRR <- (1 - RR1_matrix) > PI * (1 - RR_all_matrix)
    Method1_linearRR   <- sum(prob_matrix * condition_linearRR)

    # ----- Method 2: Pr(RR_j < 1 for all j) = Pr(hat_lambda_j < lambda0 for all j) -----
    # For each region j: Pr(Y_j <= floor(Nj * lambda0 - eps))
    # where Y_j ~ NB(mu = Nj * lambda, size = Nj * dispersion)
    Method2_probs <- numeric(J)
    for (j in seq_len(J)) {
      mu_j    <- Nj[j] * lambda
      size_j  <- Nj[j] * dispersion
      y_thr   <- floor(Nj[j] * lambda0 - 1e-10)
      Method2_probs[j] <- if (y_thr < 0) {
        0
      } else {
        stats::pnbinom(y_thr, mu = mu_j, size = size_j)
      }
    }
    Method2 <- prod(Method2_probs)

    nsim_out <- NULL

  } else {
    # ----- Monte Carlo Simulation -----
    set.seed(seed)

    # Generate total counts for all regions: Y_j ~ NB(mu = Nj*lambda, size = Nj*dispersion)
    # Dimensions: nsim rows x J columns
    yj_total_mat <- matrix(NA_integer_, nrow = nsim, ncol = J)
    for (j in seq_len(J)) {
      yj_total_mat[, j] <- stats::rnbinom(nsim,
                                          mu   = Nj[j] * lambda,
                                          size = Nj[j] * dispersion)
    }

    # Regional mean counts per patient: hat_lambda_j = Y_j / Nj[j]
    hat_lambdaj_mat <- sweep(yj_total_mat, 2, Nj, "/")

    # Overall mean count per patient: hat_lambda = sum(Y_j) / N
    hat_lambda <- rowSums(yj_total_mat) / N

    # Rate ratios
    RR1    <- hat_lambdaj_mat[, 1] / lambda0
    RR_all <- hat_lambda / lambda0

    # Method 1: Log-RR
    # Degenerate case (both numerator and overall count = 0): condition set to FALSE
    both_zero_sim  <- (yj_total_mat[, 1] == 0) & (rowSums(yj_total_mat) == 0)
    log_RR1_sim    <- ifelse(RR1    > 0, log(RR1),    -Inf)
    log_RR_all_sim <- ifelse(RR_all > 0, log(RR_all), -Inf)
    Method1_logRR  <- mean((!both_zero_sim) & (log_RR1_sim < PI * log_RR_all_sim))

    # Method 1: Linear-RR
    Method1_linearRR <- mean((1 - RR1) > PI * (1 - RR_all))

    # Method 2: Pr[hat_lambda_j < lambda0 for all j]
    Method2 <- mean(rowSums(hat_lambdaj_mat < lambda0) == J)

    nsim_out <- nsim
  }

  # ========== Output ==========
  result <- list(
    approach        = approach,
    nsim            = nsim_out,
    lambda          = lambda,
    lambda0         = lambda0,
    dispersion      = dispersion,
    Nj              = Nj,
    PI              = PI,
    Method1_logRR   = Method1_logRR,
    Method1_linearRR = Method1_linearRR,
    Method2         = Method2
  )
  class(result) <- "rcp1armCount"
  return(result)
}


#' Print Method for rcp1armCount Objects
#'
#' @param x An object of class \code{"rcp1armCount"}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
print.rcp1armCount <- function(x, ...) {
  cat("\nRegional Consistency Probability for Single-Arm MRCT\n")
  cat("Endpoint : Count (Negative Binomial)\n\n")

  if (x$approach == "simulation") {
    cat(sprintf("   Approach       : Simulation-Based (nsim = %d)\n", x$nsim))
  } else {
    cat("   Approach       : Exact Solution\n")
  }

  cat(sprintf("   Expected Count : lambda     = %.6f\n", x$lambda))
  cat(sprintf("   Control Count  : lambda0    = %.6f\n", x$lambda0))
  cat(sprintf("   Dispersion     : dispersion = %.6f\n", x$dispersion))
  cat(sprintf("   Sample Size    : Nj         = (%s)\n", paste(x$Nj, collapse = ", ")))
  cat(sprintf("   Total Size     : N          = %d\n",   sum(x$Nj)))
  cat(sprintf("   Threshold      : PI         = %.4f\n", x$PI))

  cat("\nConsistency Probabilities:\n")
  cat("   Method 1 (Region 1 vs Overall):\n")
  cat(sprintf("      Log-RR based    : %.4f\n", x$Method1_logRR))
  cat(sprintf("      Linear-RR based : %.4f\n", x$Method1_linearRR))
  cat("   Method 2 (All Regions Show Benefit):\n")
  cat(sprintf("      RR < 1          : %.4f\n", x$Method2))
  cat("\n")

  invisible(x)
}