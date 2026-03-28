#' Regional Consistency Probability for Single-Arm MRCT (Time-to-Event Endpoint)
#'
#' @description
#' Calculate the regional consistency probability (RCP) for time-to-event endpoints
#' using the hazard ratio (HR) in single-arm multi-regional clinical trials (MRCTs)
#' using the Effect Retention Approach (ERA).
#'
#' Event times are modelled by the exponential distribution with hazard rate
#' \eqn{\lambda} (treatment) relative to a known historical control hazard
#' \eqn{\lambda_0}. The treatment effect is expressed as a hazard ratio:
#' \eqn{HR = \lambda / \lambda_0 < 1} (benefit). Two effect scales are considered:
#' \itemize{
#'   \item Log-HR scale: \eqn{\log(\widehat{HR}_j) = \log(\hat{\lambda}_j / \lambda_0)}.
#'   \item Linear-HR scale: \eqn{1 - \widehat{HR}_j = 1 - \hat{\lambda}_j / \lambda_0}.
#' }
#'
#' Two evaluation methods are supported (for each scale):
#' \itemize{
#'   \item Method 1: Effect retention approach. Evaluates whether Region 1 retains
#'     at least a fraction PI of the overall treatment effect.
#'     Log-HR: \eqn{\log(\widehat{HR}_1) < \pi \times \log(\widehat{HR})};
#'     Linear-HR: \eqn{(1 - \widehat{HR}_1) > \pi \times (1 - \widehat{HR})}.
#'   \item Method 2: Simultaneous benefit across all regions.
#'     Evaluates whether all regional hazard ratios are below 1:
#'     \eqn{\widehat{HR}_j < 1} for all \eqn{j}.
#'     (Equivalent for both log-HR and linear-HR scales.)
#' }
#'
#' Two calculation approaches are available:
#' \itemize{
#'   \item \code{"formula"}: Closed-form solution based on normal approximation for
#'     \eqn{\log(\widehat{HR})} and the delta method for the linear-HR scale
#'     (Hayashi and Itoh 2018). Method 1 uses a two-block decomposition (Region 1
#'     vs regions 2..J combined), which is valid for \eqn{J \geq 2}.
#'     Method 2 supports \eqn{J \geq 2} regions.
#'   \item \code{"simulation"}: Monte Carlo simulation using individual patient data
#'     with person-years estimation of the hazard rate. Supports \eqn{J \geq 2} regions.
#' }
#'
#' @param lambda Numeric scalar. True hazard rate under the alternative hypothesis.
#'   Must be positive. Under exponential distribution, median survival =
#'   \eqn{\log(2) / \lambda}.
#' @param lambda0 Numeric scalar. Known hazard rate for the historical control
#'   (null hypothesis reference value). Must be positive.
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
#'   closed-form solution or \code{"simulation"} for Monte Carlo simulation.
#'   Default is \code{"formula"}.
#' @param nsim Positive integer. Number of Monte Carlo iterations. Used only when
#'   \code{approach = "simulation"}. Default is \code{10000}.
#' @param seed Non-negative integer. Random seed for reproducibility. Used only
#'   when \code{approach = "simulation"}. Default is \code{1}.
#'
#' @return An object of class \code{"rcp1armHazardRatio"}, which is a list containing:
#' \describe{
#'   \item{\code{approach}}{Calculation approach used (\code{"formula"} or
#'     \code{"simulation"}).}
#'   \item{\code{nsim}}{Number of Monte Carlo iterations (\code{NULL} for
#'     \code{"formula"} approach).}
#'   \item{\code{lambda}}{True hazard rate under the alternative hypothesis.}
#'   \item{\code{lambda0}}{Historical control hazard rate.}
#'   \item{\code{Nj}}{Sample sizes for each region.}
#'   \item{\code{t_a}}{Accrual period.}
#'   \item{\code{t_f}}{Follow-up period.}
#'   \item{\code{tau}}{Total study duration (\eqn{\tau = t_a + t_f}).}
#'   \item{\code{lambda_dropout}}{Dropout hazard rate (\code{NA} if \code{NULL}).}
#'   \item{\code{PI}}{Effect retention threshold.}
#'   \item{\code{Method1_logHR}}{RCP using Method 1 (log-HR scale).}
#'   \item{\code{Method1_linearHR}}{RCP using Method 1 (linear-HR scale).}
#'   \item{\code{Method2}}{RCP using Method 2 (all regions show benefit;
#'     identical for log-HR and linear-HR scales).}
#' }
#'
#' @references
#' Hayashi R, Itoh Y (2018). A reexamination of Japanese sample size calculation
#' for multiregional clinical trial evaluating survival endpoint.
#' \emph{Pharmaceutical Statistics}, 17(1): 46--55.
#'
#' Wu J (2015). Sample size calculation for the one-sample log-rank test.
#' \emph{Pharmaceutical Statistics}, 14(1): 26--33.
#'
#' @importFrom stats runif rexp pnorm
#'
#' @examples
#' # Example 1: Closed-form solution with N = 100, Region 1 has 10 subjects
#' result1 <- rcp1armHazardRatio(
#'   lambda         = log(2) / 10,
#'   lambda0        = log(2) / 5,
#'   Nj             = c(10, 90),
#'   t_a            = 3,
#'   t_f            = 10,
#'   lambda_dropout = NULL,
#'   PI             = 0.5,
#'   approach       = "formula"
#' )
#' print(result1)
#'
#' # Example 2: Monte Carlo simulation with N = 100, Region 1 has 10 subjects
#' result2 <- rcp1armHazardRatio(
#'   lambda         = log(2) / 10,
#'   lambda0        = log(2) / 5,
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
rcp1armHazardRatio <- function(lambda,
                               lambda0,
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
  if (!is.numeric(lambda0) || length(lambda0) != 1 || lambda0 <= 0) {
    stop("lambda0 must be a single positive number")
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
  N     <- sum(Nj)
  J     <- length(Nj)
  f     <- Nj[1] / N
  tau   <- t_a + t_f
  HR    <- lambda / lambda0
  delta <- log(HR)   # log-scale treatment effect; negative when lambda < lambda0

  # ========== Calculation ==========
  if (approach == "formula") {
    # Method 1: Closed-form based on normal approximation for log(HR_j).
    # Evt2 = sum(Nj[-1]) * phi represents the expected events for regions 2..J
    # combined, treating them as one block. This two-block decomposition
    # (Region 1 vs all others) is valid for J >= 2.
    # Method 2: Product of independent regional probabilities, valid for J >= 2.

    # ----- Expected event probability (Wu 2015) -----
    # phi: probability that a patient experiences an event during the study
    lambda_dropout_value <- if (is.null(lambda_dropout)) 0 else lambda_dropout
    lambda_combined      <- lambda + lambda_dropout_value

    phi <- (lambda / lambda_combined) *
      (1 - (exp(-lambda_combined * t_f) - exp(-lambda_combined * tau)) /
         (lambda_combined * t_a))

    # Expected number of events by region
    Evt1    <- Nj[1] * phi
    Evt2    <- sum(Nj[-1]) * phi   # Regions 2..J combined
    Evt_total <- N * phi

    # ----- Method 1: Log-HR -----
    # D = log(HR_1) - PI * log(HR_all), consistency requires D < 0
    # log(HR_j) ~ N(delta, 1/Evt_j) under the exponential model
    # E[D] = (1 - PI) * delta
    # Var[D] = (1 - PI*f)^2 / Evt1 + (PI*(1-f))^2 / Evt2
    mean_d  <- (1 - PI) * delta
    var_d   <- (1 - PI * f)^2 / Evt1 + (PI * (1 - f))^2 / Evt2
    sd_d    <- sqrt(var_d)

    # Pr(D < 0) = Phi(-E[D] / SD[D])
    Method1_logHR <- stats::pnorm(-mean_d / sd_d)

    # ----- Method 1: Linear-HR (Hayashi and Itoh 2018) -----
    # Consistency condition: (1 - HR_1) > PI * (1 - HR_all)
    # Transformation: g = log(HR_1) - log(1 - PI + PI * HR_all)
    # E[g]    = log(HR) - log(1 - PI + PI * HR)
    # Var[g]  = (1 / Evt_total) *
    #           [ {1 - PI + PI*(1-f)*HR}^2 / (f * {1 - PI + PI*HR}^2)
    #           + {PI*(1-f)*HR}^2 / ((1-f) * {1 - PI + PI*HR}^2) ]
    mu_g       <- log(HR) - log(1 - PI + PI * HR)
    sigma_g_sq <- (1 / Evt_total) *
      ((1 - PI + PI * (1 - f) * HR)^2 / (f       * (1 - PI + PI * HR)^2) +
       (PI * (1 - f) * HR)^2            / ((1 - f) * (1 - PI + PI * HR)^2))
    sigma_g <- sqrt(sigma_g_sq)

    # Pr(g < 0) = Phi(-E[g] / SD[g])
    Method1_linearHR <- stats::pnorm(-mu_g / sigma_g)

    # ----- Method 2: Pr(HR_j < 1 for all j) -----
    # log(HR_j) ~ N(delta, 1/Evt_j); Pr(log(HR_j) < 0) = Phi(-delta * sqrt(Evt_j))
    Method2_probs <- numeric(J)
    for (j in seq_len(J)) {
      Evt_j <- Nj[j] * phi
      Method2_probs[j] <- if (Evt_j > 0) {
        stats::pnorm(-delta * sqrt(Evt_j))
      } else {
        0
      }
    }
    Method2 <- prod(Method2_probs)

    nsim_out <- NULL

  } else {
    # ----- Monte Carlo Simulation (Individual Patient Data) -----
    set.seed(seed)

    # Matrices to accumulate events and person-years (nsim x J)
    total_events_mat <- matrix(0L,  nrow = nsim, ncol = J)
    total_py_mat     <- matrix(0.0, nrow = nsim, ncol = J)

    for (j in seq_len(J)) {
      n_total <- nsim * Nj[j]

      # Accrual times: uniform over [0, t_a]
      entry <- stats::runif(n_total, min = 0, max = t_a)

      # Event times: exponential with rate lambda
      event_time <- stats::rexp(n_total, rate = lambda)

      # Dropout times (if applicable)
      if (!is.null(lambda_dropout)) {
        dropout_time  <- stats::rexp(n_total, rate = lambda_dropout)
        true_time     <- pmin(event_time, dropout_time)
        is_event_type <- as.integer(event_time <= dropout_time)
      } else {
        true_time     <- event_time
        is_event_type <- rep(1L, n_total)
      }

      # Administrative censoring at tau
      follow_up_limit <- tau - entry
      obs_time        <- pmin(true_time, follow_up_limit)
      is_event        <- as.integer(true_time <= follow_up_limit) * is_event_type

      # Reshape to (nsim x Nj[j]) and aggregate within each simulation
      total_events_mat[, j] <- rowSums(matrix(is_event, nrow = nsim, ncol = Nj[j]))
      total_py_mat[, j]     <- rowSums(matrix(obs_time,  nrow = nsim, ncol = Nj[j]))
    }

    # ----- Hazard rate estimation: lambda_hat_j = Events_j / Person-Years_j -----
    # Person-years are always positive (observed times > 0), so no epsilon needed.
    hat_lambda_j_mat <- total_events_mat / total_py_mat
    HR_hat_j_mat     <- hat_lambda_j_mat / lambda0

    # Overall hazard: sum of all events / sum of all person-years
    hat_lambda_all <- rowSums(total_events_mat) / rowSums(total_py_mat)
    HR_hat_all     <- hat_lambda_all / lambda0

    # Method 1: Log-HR
    # When both HR_1 and HR_all are exactly 0 (zero events in all regions),
    # log(0) = -Inf on both sides: -Inf < PI * (-Inf) is NaN in R.
    # Treat this degenerate case as FALSE (insufficient data to assess consistency).
    both_zero_sim  <- (rowSums(total_events_mat) == 0)
    log_HR1_sim    <- log(HR_hat_j_mat[, 1])
    log_HR_all_sim <- log(HR_hat_all)
    Method1_logHR  <- mean((!both_zero_sim) &
                             (log_HR1_sim < PI * log_HR_all_sim))

    # Method 1: Linear-HR
    Method1_linearHR <- mean((1 - HR_hat_j_mat[, 1]) > PI * (1 - HR_hat_all))

    # Method 2: Pr[HR_j < 1 for all j]
    Method2 <- mean(rowSums(HR_hat_j_mat < 1) == J)

    nsim_out <- nsim
  }

  # ========== Output ==========
  result <- list(
    approach       = approach,
    nsim           = nsim_out,
    lambda         = lambda,
    lambda0        = lambda0,
    Nj             = Nj,
    t_a            = t_a,
    t_f            = t_f,
    tau            = tau,
    lambda_dropout = if (is.null(lambda_dropout)) NA_real_ else lambda_dropout,
    PI             = PI,
    Method1_logHR   = Method1_logHR,
    Method1_linearHR = Method1_linearHR,
    Method2         = Method2
  )
  class(result) <- "rcp1armHazardRatio"
  return(result)
}


#' Print Method for rcp1armHazardRatio Objects
#'
#' @param x An object of class \code{"rcp1armHazardRatio"}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
print.rcp1armHazardRatio <- function(x, ...) {
  cat("\nRegional Consistency Probability for Single-Arm MRCT\n")
  cat("Endpoint : Time-to-Event (Hazard Ratio)\n\n")

  if (x$approach == "simulation") {
    cat(sprintf("   Approach       : Simulation-Based (nsim = %d)\n", x$nsim))
  } else {
    cat("   Approach       : Closed-Form Solution\n")
  }

  cat(sprintf("   True Hazard    : lambda  = %.6f\n", x$lambda))
  cat(sprintf("   Control Hazard : lambda0 = %.6f\n", x$lambda0))
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

  cat("\nConsistency Probabilities:\n")
  cat("   Method 1 (Region 1 vs Overall):\n")
  cat(sprintf("      Log-HR based    : %.4f\n", x$Method1_logHR))
  cat(sprintf("      Linear-HR based : %.4f\n", x$Method1_linearHR))
  cat("   Method 2 (All Regions Show Benefit):\n")
  cat(sprintf("      HR < 1          : %.4f\n", x$Method2))
  cat("\n")

  invisible(x)
}