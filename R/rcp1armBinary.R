#' Regional Consistency Probability for Single-Arm MRCT (Binary Endpoint)
#'
#' @description
#' Calculate the regional consistency probability (RCP) for binary endpoints
#' in single-arm multi-regional clinical trials (MRCTs) using the Effect Retention
#' Approach (ERA).
#'
#' Two evaluation methods are supported:
#' \itemize{
#'   \item Method 1: Effect retention approach. Evaluates whether Region 1 retains
#'     at least a fraction PI of the overall treatment effect:
#'     \eqn{Pr[(\hat{p}_1 - p_0) > \pi \times (\hat{p} - p_0)]}.
#'   \item Method 2: Simultaneous positivity across all regions.
#'     Evaluates whether all regional response rates exceed the null value:
#'     \eqn{Pr[\hat{p}_j > p_0 \text{ for all } j]}.
#' }
#'
#' Two calculation approaches are available:
#' \itemize{
#'   \item \code{"formula"}: Exact closed-form solution via full enumeration of
#'     the binomial joint distribution. Method 1 uses a two-block decomposition
#'     (Region 1 vs regions 2..J combined), which is valid for \eqn{J \geq 2}.
#'     Method 2 supports \eqn{J \geq 2} regions.
#'   \item \code{"simulation"}: Monte Carlo simulation. Supports \eqn{J \geq 2} regions.
#' }
#'
#' @param p Numeric scalar. True response rate under the alternative hypothesis.
#'   Must be in \eqn{(0, 1)}.
#' @param p0 Numeric scalar. Null hypothesis response rate (baseline or historical
#'   control). Must be in \eqn{[0, 1)}.
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
#' @return An object of class \code{"rcp1armBinary"}, which is a list containing:
#' \describe{
#'   \item{\code{approach}}{Calculation approach used (\code{"formula"} or
#'     \code{"simulation"}).}
#'   \item{\code{nsim}}{Number of Monte Carlo iterations (\code{NULL} for
#'     \code{"formula"} approach).}
#'   \item{\code{p}}{True response rate under the alternative hypothesis.}
#'   \item{\code{p0}}{Null hypothesis response rate.}
#'   \item{\code{Nj}}{Sample sizes for each region.}
#'   \item{\code{PI}}{Effect retention threshold.}
#'   \item{\code{Method1}}{RCP using Method 1 (effect retention).}
#'   \item{\code{Method2}}{RCP using Method 2 (all regions positive).}
#' }
#'
#' @importFrom stats rbinom dbinom pbinom
#'
#' @examples
#' # Example 1: Exact solution with N = 100, Region 1 has 10 subjects
#' result1 <- rcp1armBinary(
#'   p  = 0.5,
#'   p0 = 0.2,
#'   Nj = c(10, 90),
#'   PI = 0.5,
#'   approach = "formula"
#' )
#' print(result1)
#'
#' # Example 2: Monte Carlo simulation with N = 100, Region 1 has 10 subjects
#' result2 <- rcp1armBinary(
#'   p    = 0.5,
#'   p0   = 0.2,
#'   Nj   = c(10, 90),
#'   PI   = 0.5,
#'   approach = "simulation",
#'   nsim = 10000,
#'   seed = 1
#' )
#' print(result2)
#'
#' @export
rcp1armBinary <- function(p,
                          p0,
                          Nj,
                          PI       = 0.5,
                          approach = "formula",
                          nsim     = 1e4,
                          seed     = 1) {

  # ========== Input Validation ==========
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1) {
    stop("p must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(p0) || length(p0) != 1 || p0 < 0 || p0 >= 1) {
    stop("p0 must be a single numeric value in [0, 1)")
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
    # y_others represents the pooled count for regions 2..J combined,
    # which follows Binomial(N - N1, p) by additivity. This two-block
    # decomposition (Region 1 vs all others) is valid for J >= 2.
    # Method 2: Product of independent regional probabilities, valid for J >= 2.

    # ----- Method 1: Exact Enumeration -----
    # Consistency condition: (hat_p1 - p0) > PI * (hat_p - p0)
    # Enumerate all (y1, y_others) combinations and sum joint probabilities

    # Region 1: y1 ~ Binomial(N1, p)
    y1_vals <- 0:Nj[1]
    p_y1    <- stats::dbinom(y1_vals, Nj[1], p)

    # Regions 2..J combined: y_others ~ Binomial(N - N1, p)
    n_others  <- N - Nj[1]
    yo_vals   <- 0:n_others
    p_yo      <- stats::dbinom(yo_vals, n_others, p)

    # Joint probability matrix via outer product
    prob_matrix <- outer(p_y1, p_yo)

    # Matrices of hat_p1 and hat_p for each (y1, y_others) combination
    y1_matrix  <- matrix(y1_vals, nrow = length(y1_vals), ncol = length(yo_vals))
    yo_matrix  <- matrix(yo_vals, nrow = length(y1_vals), ncol = length(yo_vals),
                         byrow = TRUE)
    hat_p1_matrix <- y1_matrix / Nj[1]
    hat_p_matrix  <- (y1_matrix + yo_matrix) / N

    # Sum probabilities satisfying the consistency condition
    condition_matrix <- (hat_p1_matrix - p0) > PI * (hat_p_matrix - p0)
    Method1 <- sum(prob_matrix * condition_matrix)

    # ----- Method 2: Product of Independent Regional Probabilities -----
    # Pr(hat_pj > p0) = Pr(Yj >= ceil(Nj * p0 + epsilon))
    # where Yj ~ Binomial(Nj, p)
    y_min_vec    <- ceiling(Nj * p0 + 1e-10)
    Method2_probs <- ifelse(
      y_min_vec > Nj,
      0,
      stats::pbinom(y_min_vec - 1, Nj, p, lower.tail = FALSE)
    )
    Method2 <- prod(Method2_probs)

    nsim_out <- NULL

  } else {
    # ----- Monte Carlo Simulation -----
    set.seed(seed)

    # Generate number of responders: Yj ~ Binomial(Nj[j], p)
    # Dimensions: nsim rows x J columns
    yj_mat <- matrix(NA_integer_, nrow = nsim, ncol = J)
    for (j in seq_len(J)) {
      yj_mat[, j] <- stats::rbinom(nsim, Nj[j], p)
    }

    # Regional response rates: hat_pj = Yj / Nj[j]
    hat_pj_mat <- sweep(yj_mat, 2, Nj, "/")

    # Overall response rate: hat_p = sum(Yj) / N
    hat_p <- rowSums(yj_mat) / N

    # Method 1: Pr[(hat_p1 - p0) > PI * (hat_p - p0)]
    Method1 <- mean((hat_pj_mat[, 1] - p0) > PI * (hat_p - p0))

    # Method 2: Pr[hat_pj > p0 for all j]
    Method2 <- mean(rowSums(hat_pj_mat > p0) == J)

    nsim_out <- nsim
  }

  # ========== Output ==========
  result <- list(
    approach = approach,
    nsim     = nsim_out,
    p        = p,
    p0       = p0,
    Nj       = Nj,
    PI       = PI,
    Method1  = Method1,
    Method2  = Method2
  )
  class(result) <- "rcp1armBinary"
  return(result)
}


#' Print Method for rcp1armBinary Objects
#'
#' @param x An object of class \code{"rcp1armBinary"}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
print.rcp1armBinary <- function(x, ...) {
  cat("\nRegional Consistency Probability for Single-Arm MRCT\n")
  cat("Endpoint : Binary\n\n")

  if (x$approach == "simulation") {
    cat(sprintf("   Approach      : Simulation-Based (nsim = %d)\n", x$nsim))
  } else {
    cat("   Approach      : Exact Solution\n")
  }

  cat(sprintf("   Response Rate : p  = %.4f\n", x$p))
  cat(sprintf("   Null Rate     : p0 = %.4f\n", x$p0))
  cat(sprintf("   Sample Size   : Nj = (%s)\n", paste(x$Nj, collapse = ", ")))
  cat(sprintf("   Total Size    : N  = %d\n",   sum(x$Nj)))
  cat(sprintf("   Threshold     : PI = %.4f\n", x$PI))

  cat("\nConsistency Probabilities:\n")
  cat(sprintf("   Method 1 (Region 1 vs Overall) : %.4f\n", x$Method1))
  cat(sprintf("   Method 2 (All Regions > p0)    : %.4f\n", x$Method2))
  cat("\n")

  invisible(x)
}