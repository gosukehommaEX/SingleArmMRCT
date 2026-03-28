#' Regional Consistency Probability for Single-Arm MRCT (RMST Endpoint)
#'
#' @description
#' Calculate the regional consistency probability (RCP) for restricted mean
#' survival time (RMST) endpoints in single-arm multi-regional clinical trials
#' (MRCTs) using the Effect Retention Approach (ERA).
#'
#' Event times are modelled by the exponential distribution with hazard rate
#' \eqn{\lambda} (treatment). The treatment effect is expressed as the difference
#' in RMST at a prespecified truncation time \eqn{\tau^*}:
#' \eqn{\delta = \mu(\tau^*) - \mu_0(\tau^*)}, where \eqn{\mu_0(\tau^*)} is
#' the historical control RMST at \eqn{\tau^*}.
#'
#' The regional RMST estimator is defined as the area under the Kaplan-Meier
#' curve up to \eqn{\tau^*}:
#' \eqn{\hat{\mu}_j(\tau^*) = \int_0^{\tau^*} \hat{S}_j(t)\,dt}.
#'
#' Two evaluation methods are supported:
#' \itemize{
#'   \item Method 1: Effect retention approach. Evaluates whether Region 1 retains
#'     at least a fraction PI of the overall treatment effect:
#'     \eqn{Pr[(\hat{\mu}_1 - \mu_0) > \pi \times (\hat{\mu} - \mu_0)]}.
#'   \item Method 2: Simultaneous benefit across all regions.
#'     Evaluates whether all regional RMST estimates exceed the historical
#'     control RMST: \eqn{Pr[\hat{\mu}_j > \mu_0 \text{ for all } j]}.
#' }
#'
#' Two calculation approaches are available:
#' \itemize{
#'   \item \code{"formula"}: Normal approximation based on exponential event times
#'     and exponential dropout. When \eqn{\tau^* \leq t_f} (the most common setting
#'     in practice), the administrative censoring survival function equals 1 on
#'     \eqn{[0, \tau^*]} and the variance integral has a closed-form solution.
#'     When \eqn{\tau^* > t_f}, the integral is evaluated numerically via
#'     \code{\link[stats]{integrate}}. Method 1 uses a two-block decomposition
#'     (Region 1 vs regions 2..J combined), valid for \eqn{J \geq 2}.
#'     Method 2 supports \eqn{J \geq 2} regions.
#'   \item \code{"simulation"}: Monte Carlo simulation. The RMST for each
#'     simulation replicate is computed as the exact area under the Kaplan-Meier
#'     step function (sum of rectangles). Supports \eqn{J \geq 2} regions.
#' }
#'
#' @details
#' \strong{Variance formula.}
#' The asymptotic variance of the RMST estimator \eqn{\hat{\mu}_j(\tau^*)} for a
#' single arm with \eqn{n} subjects is:
#' \deqn{
#'   \mathrm{Var}(\hat{\mu}_j(\tau^*)) \approx
#'   \frac{1}{n} \int_0^{\tau^*}
#'   \frac{\left(e^{-\lambda t} - e^{-\lambda \tau^*}\right)^2}
#'        {\lambda \cdot e^{-\lambda t} \cdot G(t)}
#'   \,dt,
#' }
#' where \eqn{G(t) = P(\mathrm{observed\ time} \geq t)} is the overall censoring
#' survival function:
#' \deqn{
#'   G(t) = e^{-\lambda_d t} \times G_a(t), \quad
#'   G_a(t) = \begin{cases}
#'     1              & 0 \leq t \leq t_f, \\
#'     (\tau-t)/t_a   & t_f < t \leq \tau.
#'   \end{cases}
#' }
#' \strong{Closed-form variance (\eqn{\tau^* \leq t_f}).}
#' When \eqn{\tau^* \leq t_f}, \eqn{G_a(t) = 1} on \eqn{[0, \tau^*]} and the
#' integrand reduces to
#' \eqn{e^{\lambda_d t}(1 - e^{-\lambda(\tau^*-t)})^2 / \lambda}.
#' Expanding and integrating term by term:
#' \deqn{
#'   v(\tau^*) = \frac{1}{\lambda}\Bigl[
#'     A(\lambda_d)
#'     - 2 e^{-\lambda\tau^*} A(\lambda_d + \lambda)
#'     + e^{-2\lambda\tau^*} A(\lambda_d + 2\lambda)
#'   \Bigr],
#' }
#' where \eqn{A(r) = (e^{r\tau^*} - 1)/r} for \eqn{r > 0} and
#' \eqn{A(0) = \tau^*}.
#' When there is no dropout (\eqn{\lambda_d = 0}) this further reduces to:
#' \deqn{
#'   v(\tau^*) = \tau^* - \frac{2(1-e^{-\lambda\tau^*})}{\lambda}
#'               + \frac{1-e^{-2\lambda\tau^*}}{2\lambda}.
#' }
#'
#' @param lambda Numeric scalar. True hazard rate under the alternative hypothesis.
#'   Must be positive.
#' @param tau_star Numeric scalar. Truncation time for RMST calculation. Must be
#'   positive and no greater than the total study duration \eqn{\tau = t_a + t_f}.
#' @param mu0 Numeric scalar. Historical control RMST at \code{tau_star}, i.e.,
#'   \eqn{\mu_0(\tau^*) = \int_0^{\tau^*} S_0(t)\,dt}. Must be positive and
#'   less than \code{tau_star}.
#' @param Nj Integer vector. Sample sizes for each region. For example,
#'   \code{c(10, 90)} indicates Region 1 has 10 subjects and Region 2 has 90
#'   subjects. All elements must be positive integers.
#' @param t_a Numeric scalar. Accrual period (patient enrollment duration).
#'   Must be positive.
#' @param t_f Numeric scalar. Follow-up period (additional follow-up after accrual
#'   ends). Must be positive.
#' @param lambda_dropout Numeric scalar or \code{NULL}. Dropout hazard rate.
#'   If \code{NULL} (default), no dropout is assumed. If specified, dropout times
#'   follow an exponential distribution with rate \code{lambda_dropout}.
#' @param PI Numeric scalar. Prespecified effect retention threshold for Method 1.
#'   Typically \eqn{\pi \geq 0.5}. Must be in \eqn{[0, 1]}. Default is \code{0.5}.
#' @param approach Character scalar. Calculation approach: \code{"formula"} for the
#'   normal approximation or \code{"simulation"} for Monte Carlo simulation.
#'   Default is \code{"formula"}.
#' @param nsim Positive integer. Number of Monte Carlo iterations. Used only when
#'   \code{approach = "simulation"}. Default is \code{10000}.
#' @param seed Non-negative integer. Random seed for reproducibility. Used only
#'   when \code{approach = "simulation"}. Default is \code{1}.
#'
#' @return An object of class \code{"rcp1armRMST"}, which is a list containing:
#' \describe{
#'   \item{\code{approach}}{Calculation approach used (\code{"formula"} or
#'     \code{"simulation"}).}
#'   \item{\code{formula_type}}{For \code{approach = "formula"}: either
#'     \code{"closed-form"} (when \eqn{\tau^* \leq t_f}) or
#'     \code{"numerical-integration"} (when \eqn{\tau^* > t_f}).
#'     \code{NULL} for simulation.}
#'   \item{\code{nsim}}{Number of Monte Carlo iterations (\code{NULL} for
#'     \code{"formula"} approach).}
#'   \item{\code{lambda}}{True hazard rate under the alternative hypothesis.}
#'   \item{\code{Nj}}{Sample sizes for each region.}
#'   \item{\code{t_a}}{Accrual period.}
#'   \item{\code{t_f}}{Follow-up period.}
#'   \item{\code{tau}}{Total study duration (\eqn{\tau = t_a + t_f}).}
#'   \item{\code{lambda_dropout}}{Dropout hazard rate (\code{NA} if \code{NULL}).}
#'   \item{\code{PI}}{Effect retention threshold.}
#'   \item{\code{tau_star}}{Truncation time for RMST.}
#'   \item{\code{mu0}}{Historical control RMST at \code{tau_star}.}
#'   \item{\code{mu_est}}{True RMST under the alternative:
#'     \eqn{\mu(\tau^*) = (1 - e^{-\lambda \tau^*}) / \lambda}.}
#'   \item{\code{Method1}}{RCP using Method 1 (effect retention).}
#'   \item{\code{Method2}}{RCP using Method 2 (all regions positive).}
#' }
#'
#' @references
#' Freidlin B, Korn EL (2021). Are restricted mean survival time methods
#' especially useful for noninferiority trials?
#' \emph{Clinical Trials}, 18(1): 1--9.
#'
#' Wu J (2015). Sample size calculation for the one-sample log-rank test.
#' \emph{Pharmaceutical Statistics}, 14(1): 26--33.
#'
#' @importFrom stats runif rexp pnorm integrate
#'
#' @examples
#' # Example 1: Closed-form solution (tau_star <= t_f) with N = 100
#' lam   <- log(2) / 10
#' lam0  <- log(2) / 5
#' tstar <- 8
#' mu0_val <- (1 - exp(-lam0 * tstar)) / lam0
#'
#' result1 <- rcp1armRMST(
#'   lambda         = lam,
#'   tau_star       = tstar,
#'   mu0            = mu0_val,
#'   Nj             = c(10, 90),
#'   t_a            = 3,
#'   t_f            = 10,
#'   lambda_dropout = NULL,
#'   PI             = 0.5,
#'   approach       = "formula"
#' )
#' print(result1)
#'
#' # Example 2: Monte Carlo simulation with N = 100
#' result2 <- rcp1armRMST(
#'   lambda         = lam,
#'   tau_star       = tstar,
#'   mu0            = mu0_val,
#'   Nj             = c(10, 90),
#'   t_a            = 3,
#'   t_f            = 10,
#'   lambda_dropout = NULL,
#'   PI             = 0.5,
#'   approach       = "simulation",
#'   nsim           = 10000,
#'   seed           = 1
#' )
#' print(result2)
#'
#' @export
rcp1armRMST <- function(lambda,
                        tau_star,
                        mu0,
                        Nj,
                        t_a,
                        t_f,
                        lambda_dropout = NULL,
                        PI             = 0.5,
                        approach       = "formula",
                        nsim           = 1e4,
                        seed           = 1) {

  # ========== Input Validation ==========
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) {
    stop("lambda must be a single positive number")
  }
  if (!is.numeric(tau_star) || length(tau_star) != 1 || tau_star <= 0) {
    stop("tau_star must be a single positive number")
  }
  if (!is.numeric(mu0) || length(mu0) != 1 || mu0 <= 0) {
    stop("mu0 must be a single positive number")
  }
  if (!is.numeric(Nj) || any(Nj <= 0) || any(Nj != as.integer(Nj))) {
    stop("Nj must be a vector of positive integers")
  }
  if (approach == "simulation" && any(Nj < 3)) {
    stop("Nj must have all elements >= 3 for the simulation approach (required by the KM RMST estimator)")
  }
  if (!is.numeric(t_a) || length(t_a) != 1 || t_a <= 0) {
    stop("t_a must be a single positive number")
  }
  if (!is.numeric(t_f) || length(t_f) != 1 || t_f <= 0) {
    stop("t_f must be a single positive number")
  }
  if (!is.null(lambda_dropout)) {
    if (!is.numeric(lambda_dropout) || length(lambda_dropout) != 1 ||
        lambda_dropout <= 0) {
      stop("lambda_dropout must be a single positive number or NULL")
    }
  }
  tau <- t_a + t_f
  if (tau_star > tau) {
    stop("tau_star must not exceed the total study duration tau = t_a + t_f")
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
  N        <- sum(Nj)
  J        <- length(Nj)
  f        <- Nj[1] / N
  lambda_d <- if (is.null(lambda_dropout)) 0 else lambda_dropout

  # True RMST under the alternative: integral of S(t) = exp(-lambda*t) from 0 to tau_star
  mu_est <- (1 - exp(-lambda * tau_star)) / lambda

  # ========== Calculation ==========
  if (approach == "formula") {

    # ----- Variance kernel v(tau_star) -----
    # Var(mu_hat_j) = v(tau_star) / n_j
    #
    # v(tau_star) = integral_0^{tau_star}
    #   (S(t) - S(tau_star))^2 / (lambda * S(t) * G(t)) dt
    # = integral_0^{tau_star}
    #   exp(lambda_d * t) * (1 - exp(-lambda*(tau_star-t)))^2 / (lambda * G_a(t)) dt
    #
    # where G_a(t) = 1 for t <= t_f, and (tau-t)/t_a for t_f < t <= tau.

    if (tau_star <= t_f) {
      # --- Closed-form solution: G_a(t) = 1 throughout [0, tau_star] ---
      #
      # Expand (1 - exp(-lambda*(tau_star-t)))^2:
      #   = 1 - 2*exp(-lambda*(tau_star-t)) + exp(-2*lambda*(tau_star-t))
      #
      # Each term integrates as exp(lambda_d * t) * exp(c * t) = exp((lambda_d+c)*t),
      # so the integral from 0 to tau_star is A(lambda_d + c) where:
      #   A(r) = (exp(r * tau_star) - 1) / r  for r > 0
      #   A(0) = tau_star
      A_fn <- function(r) {
        if (abs(r) < .Machine$double.eps^0.5) tau_star else (exp(r * tau_star) - 1) / r
      }

      v_rmst <- (1 / lambda) * (
        A_fn(lambda_d) -
          2 * exp(-lambda  * tau_star) * A_fn(lambda_d + lambda) +
          exp(-2 * lambda * tau_star) * A_fn(lambda_d + 2 * lambda)
      )
      formula_type <- "closed-form"

    } else {
      # --- Numerical integration: tau_star > t_f ---
      # Split at t_f to handle the kink in G_a(t).
      integrand <- function(t) {
        G_a         <- ifelse(t <= t_f, 1, (tau - t) / t_a)
        numerator   <- exp(lambda_d * t) * (1 - exp(-lambda * (tau_star - t)))^2
        denominator <- lambda * G_a
        ifelse(G_a > 0, numerator / denominator, 0)
      }

      v1 <- stats::integrate(integrand, lower = 0,   upper = t_f,
                             subdivisions = 500L,
                             rel.tol = .Machine$double.eps^0.5)$value
      v2 <- stats::integrate(integrand, lower = t_f, upper = tau_star,
                             subdivisions = 500L,
                             rel.tol = .Machine$double.eps^0.5)$value
      v_rmst <- v1 + v2
      formula_type <- "numerical-integration"
    }

    # Variance per region
    var_mu_j       <- v_rmst / Nj          # vector of length J
    var_mu_region1 <- var_mu_j[1]
    var_mu_region2 <- v_rmst / sum(Nj[-1]) # regions 2..J treated as one block

    # ----- Method 1 -----
    # D = (hat_mu_1 - mu0) - PI * (hat_mu_all - mu0)
    # Under homogeneity: E[D] = (1 - PI) * (mu_est - mu0)
    # Var[D] = (1 - PI*f)^2 * Var(mu_hat_1) + (PI*(1-f))^2 * Var(mu_hat_2)
    mean_d <- (1 - PI) * (mu_est - mu0)
    var_d  <- (1 - PI * f)^2 * var_mu_region1 +
      (PI * (1 - f))^2 * var_mu_region2
    sd_d   <- sqrt(var_d)

    Method1 <- stats::pnorm(mean_d / sd_d)

    # ----- Method 2 -----
    # Pr(hat_mu_j > mu0 for all j) = product_j Phi((mu_est - mu0) / sqrt(Var_j))
    Method2_probs <- stats::pnorm((mu_est - mu0) / sqrt(var_mu_j))
    Method2       <- prod(Method2_probs)

    nsim_out <- NULL

  } else {
    # ----- Monte Carlo Simulation -----
    set.seed(seed)

    # Exact RMST via KM step-function rectangle sum up to tau_star.
    #
    # KM is a left-continuous step function dropping only at event times.
    # The exact area from 0 to tau_star is:
    #   RMST = sum_i S(t_i-) * (min(t_i, tau_star) - min(t_{i-1}, tau_star))
    # where t_i are all ordered observed times (events and censored), t_0 = 0,
    # and S(t_i-) is the KM survival just before the i-th observation.
    # Censored observations contribute step factor 1 (no drop in S) but are
    # retained to maintain correct risk-set sizes.
    #
    # Vectorization: flatten row-major, sort globally by (row_id, time), reshape.
    # Log-cumulative-sum of KM step factors uses row-wise apply(cumsum) on the
    # interior columns, then prepends a zero column for S(t_1-) = 1.
    fast_rmst <- function(X_mat, D_mat, tau_star) {
      nsim_local <- nrow(X_mat)
      n_pts      <- ncol(X_mat)

      # Flatten and sort by (row index, observed time)
      X_flat     <- as.vector(t(X_mat))
      D_flat     <- as.vector(t(D_mat))
      row_id     <- rep(seq_len(nsim_local), each = n_pts)
      global_ord <- order(row_id, X_flat)

      X_sorted <- matrix(X_flat[global_ord], nrow = nsim_local, byrow = TRUE)
      D_sorted <- matrix(D_flat[global_ord], nrow = nsim_local, byrow = TRUE)

      # Risk set sizes: n_pts, n_pts-1, ..., 1 (identical across all rows)
      n_risk_mat <- matrix(n_pts:1, nrow = nsim_local, ncol = n_pts, byrow = TRUE)

      # KM log-step factors: log(1 - D_i / n_i)
      # Events: log(1 - 1/n_i) < 0; censored: log(1) = 0 (no drop)
      log_step_mat <- log(1 - D_sorted / n_risk_mat)

      # S(t_i-): KM survival just BEFORE the i-th ordered observation.
      # S(t_1-) = 1 (prepend 0 in log scale).
      # S(t_i-) = exp(cumsum of log_step up to column i-1).
      cs_mat        <- t(apply(log_step_mat[, -n_pts, drop = FALSE], 1, cumsum))
      log_S_pre_mat <- cbind(0, cs_mat)
      S_pre_mat     <- exp(log_S_pre_mat)   # nsim x n_pts

      # Rectangle widths: dt_i = min(t_i, tau_star) - min(t_{i-1}, tau_star) >= 0
      t_trunc <- pmin(X_sorted, tau_star)
      t_prev  <- cbind(0, t_trunc[, -n_pts, drop = FALSE])
      dt_mat  <- t_trunc - t_prev

      # Exact RMST = sum_i S(t_i-) * dt_i
      rowSums(S_pre_mat * dt_mat)
    }

    # Generate patient data and RMST estimates for each region
    mu_hat_j_mat <- matrix(NA_real_, nrow = nsim, ncol = J)
    X_all_list   <- vector("list", J)
    D_all_list   <- vector("list", J)

    for (j in seq_len(J)) {
      # Accrual times: uniform over [0, t_a]
      A_mat <- matrix(stats::runif(nsim * Nj[j], min = 0, max = t_a),
                      nrow = nsim, ncol = Nj[j])

      # Event times: exponential with rate lambda
      T_true_mat <- matrix(stats::rexp(nsim * Nj[j], rate = lambda),
                           nrow = nsim, ncol = Nj[j])

      # Administrative censoring time for each patient
      C_mat <- tau - A_mat

      # Dropout censoring (if applicable)
      if (!is.null(lambda_dropout)) {
        C_dropout_mat <- matrix(stats::rexp(nsim * Nj[j], rate = lambda_dropout),
                                nrow = nsim, ncol = Nj[j])
        C_mat <- pmin(C_mat, C_dropout_mat)
      }

      # Observed times and event indicators
      X_mat <- pmin(T_true_mat, C_mat)
      D_mat <- (T_true_mat <= C_mat) * 1L

      mu_hat_j_mat[, j] <- fast_rmst(X_mat, D_mat, tau_star)
      X_all_list[[j]]   <- X_mat
      D_all_list[[j]]   <- D_mat
    }

    # Overall RMST estimate (all regions combined)
    X_all_mat  <- do.call(cbind, X_all_list)
    D_all_mat  <- do.call(cbind, D_all_list)
    mu_hat_all <- fast_rmst(X_all_mat, D_all_mat, tau_star)

    # Method 1: Pr[(hat_mu_1 - mu0) > PI * (hat_mu_all - mu0)]
    Method1 <- mean((mu_hat_j_mat[, 1] - mu0) > PI * (mu_hat_all - mu0))

    # Method 2: Pr[hat_mu_j > mu0 for all j]
    Method2 <- mean(rowSums(mu_hat_j_mat > mu0) == J)

    nsim_out     <- nsim
    formula_type <- NULL
  }

  # ========== Output ==========
  result <- list(
    approach       = approach,
    formula_type   = formula_type,
    nsim           = nsim_out,
    lambda         = lambda,
    Nj             = Nj,
    t_a            = t_a,
    t_f            = t_f,
    tau            = tau,
    lambda_dropout = if (is.null(lambda_dropout)) NA_real_ else lambda_dropout,
    PI             = PI,
    tau_star       = tau_star,
    mu0            = mu0,
    mu_est         = mu_est,
    Method1        = Method1,
    Method2        = Method2
  )
  class(result) <- "rcp1armRMST"
  return(result)
}


#' Print Method for rcp1armRMST Objects
#'
#' @param x An object of class \code{"rcp1armRMST"}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
print.rcp1armRMST <- function(x, ...) {
  cat("\nRegional Consistency Probability for Single-Arm MRCT\n")
  cat("Endpoint : Restricted Mean Survival Time (RMST)\n\n")

  if (x$approach == "simulation") {
    cat(sprintf("   Approach       : Simulation-Based (nsim = %d)\n", x$nsim))
  } else {
    label <- if (identical(x$formula_type, "closed-form")) {
      "Closed-Form Solution"
    } else {
      "Normal Approximation (Numerical Integration)"
    }
    cat(sprintf("   Approach       : %s\n", label))
  }

  cat(sprintf("   True Hazard    : lambda  = %.6f\n", x$lambda))
  cat(sprintf("   Sample Size    : Nj      = (%s)\n", paste(x$Nj, collapse = ", ")))
  cat(sprintf("   Total Size     : N       = %d\n",   sum(x$Nj)))
  cat(sprintf("   Accrual Period : t_a     = %.2f\n", x$t_a))
  cat(sprintf("   Follow-up      : t_f     = %.2f\n", x$t_f))
  cat(sprintf("   Study Duration : tau     = %.2f\n", x$tau))
  if (is.na(x$lambda_dropout)) {
    cat("   Dropout Hazard : lambda_d = NA\n")
  } else {
    cat(sprintf("   Dropout Hazard : lambda_d = %.6f\n", x$lambda_dropout))
  }
  cat(sprintf("   Threshold      : PI      = %.4f\n", x$PI))
  cat(sprintf("   Trunc. Time    : tau*    = %.2f\n", x$tau_star))
  cat(sprintf("   Control RMST   : mu0     = %.4f\n", x$mu0))
  cat(sprintf("   True RMST      : mu_est  = %.4f\n", x$mu_est))

  cat("\nConsistency Probabilities:\n")
  cat(sprintf("   Method 1 (Region 1 vs Overall)  : %.4f\n", x$Method1))
  cat(sprintf("   Method 2 (All Regions > mu0)    : %.4f\n", x$Method2))
  cat("\n")

  invisible(x)
}
