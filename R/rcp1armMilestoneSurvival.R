#' Regional Consistency Probability for Single-Arm MRCT (Milestone Survival Endpoint)
#'
#' @description
#' Calculate the regional consistency probability (RCP) for milestone survival
#' endpoints in single-arm multi-regional clinical trials (MRCTs) using the Effect
#' Retention Approach (ERA).
#'
#' Event times are modelled by the exponential distribution with hazard rate
#' \eqn{\lambda} (treatment). The treatment effect at a prespecified evaluation
#' time \eqn{t_{\mathrm{eval}}} is expressed as the difference in survival rates:
#' \eqn{\delta = S(t_{\mathrm{eval}}) - S_0(t_{\mathrm{eval}})}, where
#' \eqn{S_0} is the historical control survival rate at \eqn{t_{\mathrm{eval}}}.
#'
#' Two evaluation methods are supported:
#' \itemize{
#'   \item Method 1: Effect retention approach. Evaluates whether Region 1 retains
#'     at least a fraction PI of the overall treatment effect:
#'     \eqn{Pr[(\hat{S}_1(t) - S_0(t)) > \pi \times (\hat{S}(t) - S_0(t))]}.
#'   \item Method 2: Simultaneous benefit across all regions.
#'     Evaluates whether all regional Kaplan-Meier estimates exceed the historical
#'     control: \eqn{Pr[\hat{S}_j(t) > S_0(t) \text{ for all } j]}.
#' }
#'
#' Two calculation approaches are available:
#' \itemize{
#'   \item \code{"formula"}: Closed-form or semi-analytical solution based on the
#'     asymptotic variance of the Kaplan-Meier estimator derived from Greenwood's
#'     formula (Tang 2022). When \eqn{t_{\mathrm{eval}} \leq t_f}, the
#'     administrative censoring survival function \eqn{G_a(t) = 1} and the
#'     variance integral has a closed-form solution. When
#'     \eqn{t_{\mathrm{eval}} > t_f}, the integral is evaluated numerically via
#'     \code{\link[stats]{integrate}}. Method 1 uses a two-block decomposition
#'     (Region 1 vs regions 2..J combined). Method 2 supports \eqn{J \geq 2}
#'     regions.
#'   \item \code{"simulation"}: Monte Carlo simulation with vectorized Kaplan-Meier
#'     estimation. Supports \eqn{J \geq 2} regions.
#' }
#'
#' @details
#' \strong{Variance formula (Tang 2022).}
#' The asymptotic variance of the Kaplan-Meier estimator \eqn{\hat{S}_j(t)} for a
#' single arm with \eqn{N_j} subjects is derived from Greenwood's formula:
#' \deqn{
#'   \mathrm{Var}[\hat{S}_j(t)] \approx
#'   \frac{S^2(t)}{N_j} \int_0^t \frac{\lambda(u)}{S(u) \cdot G(u)} \, du,
#' }
#' where \eqn{G(u) = e^{-\lambda_d u} \cdot G_a(u)} is the overall censoring
#' survival function, and
#' \deqn{
#'   G_a(u) = \begin{cases}
#'     1              & 0 \leq u \leq t_f, \\
#'     (\tau - u)/t_a & t_f < u \leq \tau.
#'   \end{cases}
#' }
#' Under the exponential model \eqn{S(u) = e^{-\lambda u}}, this simplifies to:
#' \deqn{
#'   \mathrm{Var}[\hat{S}_j(t)] \approx
#'   \frac{e^{-2\lambda t}}{N_j} \int_0^t
#'   \frac{\lambda \, e^{(\lambda + \lambda_d) u}}{G_a(u)} \, du.
#' }
#'
#' \strong{Closed-form solution (\eqn{t_{\mathrm{eval}} \leq t_f}).}
#' When \eqn{t_{\mathrm{eval}} \leq t_f}, \eqn{G_a(u) = 1} throughout
#' \eqn{[0, t_{\mathrm{eval}}]}, and the integral reduces to:
#' \deqn{
#'   \int_0^t \lambda \, e^{(\lambda + \lambda_d) u} \, du
#'   = \frac{\lambda}{\lambda + \lambda_d}
#'     \bigl( e^{(\lambda + \lambda_d) t} - 1 \bigr).
#' }
#' Therefore:
#' \deqn{
#'   \mathrm{Var}[\hat{S}_j(t)] =
#'   \frac{e^{-2\lambda t}}{N_j} \cdot
#'   \frac{\lambda}{\lambda + \lambda_d}
#'   \bigl( e^{(\lambda + \lambda_d) t} - 1 \bigr).
#' }
#' When there is no dropout (\eqn{\lambda_d = 0}), this further simplifies to the
#' binomial variance:
#' \deqn{
#'   \mathrm{Var}[\hat{S}_j(t)] = \frac{S(t)(1 - S(t))}{N_j}.
#' }
#'
#' @param lambda Numeric scalar. True hazard rate under the alternative hypothesis.
#'   Must be positive.
#' @param t_eval Numeric scalar. Evaluation time point for the milestone survival
#'   probability. Must be positive.
#' @param S0 Numeric scalar. Historical control survival rate at \code{t_eval}.
#'   Must be in \eqn{(0, 1]}.
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
#'   closed-form (or semi-analytical) solution or \code{"simulation"} for Monte
#'   Carlo simulation. Default is \code{"formula"}.
#' @param nsim Positive integer. Number of Monte Carlo iterations. Used only when
#'   \code{approach = "simulation"}. Default is \code{10000}.
#' @param seed Non-negative integer. Random seed for reproducibility. Used only
#'   when \code{approach = "simulation"}. Default is \code{1}.
#'
#' @return An object of class \code{"rcp1armMilestoneSurvival"}, which is a list
#' containing:
#' \describe{
#'   \item{\code{approach}}{Calculation approach used (\code{"formula"} or
#'     \code{"simulation"}).}
#'   \item{\code{formula_type}}{For \code{approach = "formula"}: either
#'     \code{"closed-form"} (when \eqn{t_{\mathrm{eval}} \leq t_f}) or
#'     \code{"numerical-integration"} (when \eqn{t_{\mathrm{eval}} > t_f}).
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
#'   \item{\code{eval_time}}{Milestone evaluation time point.}
#'   \item{\code{S0}}{Historical control survival rate at \code{eval_time}.}
#'   \item{\code{S_est}}{True survival rate under the alternative at
#'     \code{eval_time}: \eqn{S(t) = e^{-\lambda t}}.}
#'   \item{\code{Method1}}{RCP using Method 1 (effect retention).}
#'   \item{\code{Method2}}{RCP using Method 2 (all regions positive).}
#' }
#'
#' @references
#' Tang Y (2022). Complex survival trial design by the product integration method.
#' \emph{Statistics in Medicine}, 41(4): 798--814.
#'
#' Wu J (2015). Sample size calculation for the one-sample log-rank test.
#' \emph{Pharmaceutical Statistics}, 14(1): 26--33.
#'
#' @importFrom stats pnorm integrate runif rexp
#'
#' @examples
#' # Example 1: Closed-form solution (t_eval <= t_f) with N = 100
#' result1 <- rcp1armMilestoneSurvival(
#'   lambda         = log(2) / 10,
#'   t_eval         = 8,
#'   S0             = exp(-log(2) * 8 / 5),
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
#' result2 <- rcp1armMilestoneSurvival(
#'   lambda         = log(2) / 10,
#'   t_eval         = 8,
#'   S0             = exp(-log(2) * 8 / 5),
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
rcp1armMilestoneSurvival <- function(lambda,
                                     t_eval,
                                     S0,
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
  if (!is.numeric(t_eval) || length(t_eval) != 1 || t_eval <= 0) {
    stop("t_eval must be a single positive number")
  }
  if (!is.numeric(S0) || length(S0) != 1 || S0 <= 0 || S0 > 1) {
    stop("S0 must be a single numeric value in (0, 1]")
  }
  if (!is.numeric(Nj) || any(Nj <= 0) || any(Nj != as.integer(Nj))) {
    stop("Nj must be a vector of positive integers")
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
  tau      <- t_a + t_f
  lambda_d <- if (is.null(lambda_dropout)) 0 else lambda_dropout
  
  # True survival rate at t_eval under the exponential model
  S_est <- exp(-lambda * t_eval)
  
  # ========== Calculation ==========
  if (approach == "formula") {
    
    # ---- Variance of the KM estimator via Greenwood's formula (Tang 2022) ----
    #
    # Var[S_hat_j(t)] = S^2(t) / N_j * integral_0^t  lambda(u) / (S(u) * G(u)) du
    #
    # Under exponential model: S(u) = exp(-lambda * u), lambda(u) = lambda,
    #   G(u) = exp(-lambda_d * u) * G_a(u),
    #   G_a(u) = 1              for u in [0, t_f]
    #           (tau - u) / t_a for u in (t_f, tau]
    #
    # This simplifies to:
    #   Var[S_hat_j(t)] = exp(-2*lambda*t) / N_j
    #                     * integral_0^t  lambda * exp((lambda + lambda_d)*u) / G_a(u) du
    #
    # Let I(t) = integral_0^t  lambda * exp((lambda + lambda_d)*u) / G_a(u) du.
    # Then Var[S_hat_j(t)] = exp(-2*lambda*t) * I(t) / N_j.
    
    r <- lambda + lambda_d   # combined exponent rate
    
    if (t_eval <= t_f) {
      # --- Closed-form: G_a(u) = 1 throughout [0, t_eval] ---
      # I(t) = lambda / r * (exp(r * t) - 1)
      #   (limit as r -> 0: I(t) = lambda * t)
      I_val <- if (abs(r) < .Machine$double.eps^0.5) {
        lambda * t_eval
      } else {
        (lambda / r) * (exp(r * t_eval) - 1)
      }
      formula_type <- "closed-form"
      
    } else {
      # --- Numerical integration: t_eval > t_f ---
      # Split at t_f where G_a(u) has a kink.
      #
      # Segment 1: u in [0, t_f]  -> G_a(u) = 1
      # I1 = lambda / r * (exp(r * t_f) - 1)
      I1 <- if (abs(r) < .Machine$double.eps^0.5) {
        lambda * t_f
      } else {
        (lambda / r) * (exp(r * t_f) - 1)
      }
      
      # Segment 2: u in (t_f, t_eval] -> G_a(u) = (tau - u) / t_a
      integrand2 <- function(u) {
        G_a <- (tau - u) / t_a
        lambda * exp(r * u) / G_a
      }
      I2 <- stats::integrate(integrand2, lower = t_f, upper = t_eval,
                             subdivisions = 500L,
                             rel.tol = .Machine$double.eps^0.5)$value
      
      I_val        <- I1 + I2
      formula_type <- "numerical-integration"
    }
    
    # Per-subject variance kernel (multiplied by N_j to get per-region variance)
    v_km <- exp(-2 * lambda * t_eval) * I_val
    
    # Per-region variance: Var[S_hat_j] = v_km / N_j
    var_S_j       <- v_km / Nj
    var_S_region1 <- var_S_j[1]
    # Regions 2..J combined: treat as one block with N_rest = sum(Nj[-1]) subjects
    var_S_region2 <- v_km / sum(Nj[-1])
    
    # ----- Method 1 -----
    # D = (S_hat_1 - S0) - PI * (S_hat_all - S0)
    #   = (1 - PI*f) * S_hat_1 - PI*(1-f) * S_hat_2 - (1 - PI) * S0
    # Under homogeneity: E[D] = (1 - PI) * (S_est - S0)
    # Var[D] = (1 - PI*f)^2 * Var[S_hat_1] + (PI*(1-f))^2 * Var[S_hat_2]
    mean_d <- (1 - PI) * (S_est - S0)
    var_d  <- (1 - PI * f)^2 * var_S_region1 +
      (PI * (1 - f))^2 * var_S_region2
    sd_d   <- sqrt(var_d)
    
    Method1 <- stats::pnorm(mean_d / sd_d)
    
    # ----- Method 2 -----
    # Pr(S_hat_j > S0 for all j) = prod_j Phi((S_est - S0) / sqrt(Var[S_hat_j]))
    Method2_probs <- stats::pnorm((S_est - S0) / sqrt(var_S_j))
    Method2       <- prod(Method2_probs)
    
    nsim_out <- NULL
    
  } else {
    # ----- Monte Carlo Simulation (Vectorized Kaplan-Meier) -----
    set.seed(seed)
    
    # Fast KM estimator at t_eval using fully vectorized operations.
    # Arguments:
    #   X_mat : nsim x n_pts matrix of observed times
    #   D_mat : nsim x n_pts matrix of event indicators (1 = event, 0 = censored)
    # Returns: numeric vector of length nsim with S_hat(t_eval) per simulation.
    #
    # Key optimizations:
    #   1. Row-wise sorting replaced by a single global order() call on the
    #      row-major flattened vector, keyed by (row index, observed time).
    #      This avoids the slow t(apply(X_mat, 1, order)) pattern.
    #   2. Row-wise product of KM step factors computed as
    #      exp(rowSums(log(prob_step))), using C-level rowSums() instead of
    #      apply(prob_step, 1, prod).
    fast_km_at_t <- function(X_mat, D_mat, t_eval) {
      nsim_local <- nrow(X_mat)
      n_pts      <- ncol(X_mat)
      
      # Flatten matrices in row-major order and sort by (row index, observed time)
      X_flat     <- as.vector(t(X_mat))
      D_flat     <- as.vector(t(D_mat))
      row_id     <- rep(seq_len(nsim_local), each = n_pts)
      global_ord <- order(row_id, X_flat)
      
      X_sorted   <- matrix(X_flat[global_ord], nrow = nsim_local, byrow = TRUE)
      D_sorted   <- matrix(D_flat[global_ord], nrow = nsim_local, byrow = TRUE)
      
      # Risk set sizes: n_pts, n_pts-1, ..., 1 (identical across all rows)
      n_risk_mat <- matrix(n_pts:1, nrow = nsim_local, ncol = n_pts, byrow = TRUE)
      
      # KM step factors: (1 - d_i / n_i) for events at or before t_eval
      at_risk_before_t <- X_sorted <= t_eval
      prob_step        <- 1 - (D_sorted * at_risk_before_t) / n_risk_mat
      
      # Product of step factors via log-sum for speed (C-level rowSums)
      exp(rowSums(log(prob_step)))
    }
    
    # Generate patient data and KM estimates for each region
    S_hat_j_mat <- matrix(NA_real_, nrow = nsim, ncol = J)
    X_all_list  <- vector("list", J)
    D_all_list  <- vector("list", J)
    
    for (j in seq_len(J)) {
      # Accrual times: uniform over [0, t_a]
      A_mat <- matrix(stats::runif(nsim * Nj[j], min = 0, max = t_a),
                      nrow = nsim, ncol = Nj[j])
      # Event times: exponential with rate lambda
      T_true_mat <- matrix(stats::rexp(nsim * Nj[j], rate = lambda),
                           nrow = nsim, ncol = Nj[j])
      # Administrative censoring
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
      
      S_hat_j_mat[, j] <- fast_km_at_t(X_mat, D_mat, t_eval)
      X_all_list[[j]]  <- X_mat
      D_all_list[[j]]  <- D_mat
    }
    
    # Overall KM estimate (all regions combined)
    X_all_mat <- do.call(cbind, X_all_list)
    D_all_mat <- do.call(cbind, D_all_list)
    S_hat_all <- fast_km_at_t(X_all_mat, D_all_mat, t_eval)
    
    # Method 1: Pr[(S_hat_1 - S0) > PI * (S_hat_all - S0)]
    Method1 <- mean((S_hat_j_mat[, 1] - S0) > PI * (S_hat_all - S0))
    
    # Method 2: Pr[S_hat_j > S0 for all j]
    Method2 <- mean(rowSums(S_hat_j_mat > S0) == J)
    
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
    eval_time      = t_eval,
    S0             = S0,
    S_est          = S_est,
    Method1        = Method1,
    Method2        = Method2
  )
  class(result) <- "rcp1armMilestoneSurvival"
  return(result)
}


#' Print Method for rcp1armMilestoneSurvival Objects
#'
#' @param x An object of class \code{"rcp1armMilestoneSurvival"}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
print.rcp1armMilestoneSurvival <- function(x, ...) {
  cat("\nRegional Consistency Probability for Single-Arm MRCT\n")
  cat("Endpoint : Milestone Survival\n\n")
  
  if (x$approach == "simulation") {
    cat(sprintf("   Approach       : Simulation-Based (nsim = %d)\n", x$nsim))
  } else {
    label <- if (identical(x$formula_type, "closed-form")) {
      "Closed-Form Solution (Greenwood)"
    } else {
      "Semi-Analytical (Greenwood + Numerical Integration)"
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
  cat(sprintf("   Eval Time      : t_eval  = %.2f\n", x$eval_time))
  cat(sprintf("   Control Surv   : S0      = %.4f\n", x$S0))
  cat(sprintf("   True Surv      : S_est   = %.4f\n", x$S_est))
  
  cat("\nConsistency Probabilities:\n")
  cat(sprintf("   Method 1 (Region 1 vs Overall)  : %.4f\n", x$Method1))
  cat(sprintf("   Method 2 (All Regions > S0)     : %.4f\n", x$Method2))
  cat("\n")
  
  invisible(x)
}