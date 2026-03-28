#' Plot Regional Consistency Probability for Single-Arm MRCT (Milestone Survival Endpoint)
#'
#' @description
#' Generate a faceted plot of Regional Consistency Probability (RCP) as a function
#' of the regional allocation proportion \eqn{f_1} for milestone survival endpoints.
#' Formula and simulation results are shown together for both Method 1 and Method 2.
#' Facet columns correspond to total sample sizes specified in \code{N_vec}.
#'
#' Regional sample sizes are allocated as:
#' \eqn{N_{j1} = \lfloor N \times f_1 \rfloor} and
#' \eqn{N_{j2} = \cdots = N_{jJ} = (N - N_{j1}) / (J - 1)}.
#'
#' @param lambda Numeric scalar. True hazard rate under the alternative hypothesis.
#'   Must be positive. Default is \code{log(2) / 10}.
#' @param t_eval Numeric scalar. Milestone evaluation time point. Must be positive.
#'   Default is \code{8}.
#' @param S0 Numeric scalar. Historical control survival rate at \code{t_eval}.
#'   Must be in \eqn{(0, 1)}. Default is \code{exp(-log(2) * 8 / 5)}.
#' @param t_a Numeric scalar. Accrual period. Must be positive. Default is \code{3}.
#' @param t_f Numeric scalar. Follow-up period. Must be positive. Default is \code{10}.
#' @param lambda_dropout Numeric scalar or \code{NULL}. Dropout hazard rate.
#'   If \code{NULL} (default), no dropout is assumed.
#' @param PI Numeric scalar. Effect retention threshold for Method 1.
#'   Must be in \eqn{[0, 1]}. Default is \code{0.5}.
#' @param N_vec Integer vector. Total sample sizes for each facet column.
#'   Default is \code{c(20, 40, 100)}.
#' @param J Positive integer (>= 2). Number of regions. Default is \code{3}.
#' @param f1_seq Numeric vector. Sequence of Region 1 allocation proportions.
#'   Each value must be in \eqn{(0, 1)}. Default is \code{seq(0.1, 0.9, by = 0.1)}.
#' @param nsim Positive integer. Number of Monte Carlo iterations for simulation.
#'   Default is \code{10000}.
#' @param seed Non-negative integer. Random seed for simulation. Default is \code{1}.
#' @param base_size Positive numeric. Base font size in points passed to
#'   \code{\link[ggplot2]{theme}}. Use larger values (e.g., \code{28}) for
#'   presentation slides and smaller values (e.g., \code{11}) for vignettes or
#'   reports. Default is \code{28}.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_hline scale_color_manual
#'   scale_linetype_manual scale_x_continuous scale_y_continuous facet_wrap
#'   labs theme_bw theme element_text element_blank unit guide_legend guides
#'   label_parsed
#'
#' @examples
#' p <- plot_rcp1armMilestoneSurvival(
#'   lambda = log(2) / 10,
#'   t_eval = 8,
#'   S0     = exp(-log(2) * 8 / 5),
#'   t_a    = 3,
#'   t_f    = 10,
#'   PI     = 0.5,
#'   N_vec  = c(20, 40, 100),
#'   J      = 3
#' )
#' print(p)
#'
#' @export
plot_rcp1armMilestoneSurvival <- function(lambda         = log(2) / 10,
                                          t_eval         = 8,
                                          S0             = exp(-log(2) * 8 / 5),
                                          t_a            = 3,
                                          t_f            = 10,
                                          lambda_dropout = NULL,
                                          PI             = 0.5,
                                          N_vec          = c(20, 40, 100),
                                          J              = 3,
                                          f1_seq         = seq(0.1, 0.9, by = 0.1),
                                          nsim           = 1e4,
                                          seed           = 1,
                                          base_size      = 28) {

  # ========== Input Validation ==========
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) {
    stop("lambda must be a single positive number")
  }
  if (!is.numeric(t_eval) || length(t_eval) != 1 || t_eval <= 0) {
    stop("t_eval must be a single positive number")
  }
  if (!is.numeric(S0) || length(S0) != 1 || S0 <= 0 || S0 >= 1) {
    stop("S0 must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(t_a) || length(t_a) != 1 || t_a <= 0) {
    stop("t_a must be a single positive number")
  }
  if (!is.numeric(t_f) || length(t_f) != 1 || t_f <= 0) {
    stop("t_f must be a single positive number")
  }
  if (!is.null(lambda_dropout) &&
      (!is.numeric(lambda_dropout) || length(lambda_dropout) != 1 ||
       lambda_dropout <= 0)) {
    stop("lambda_dropout must be a single positive number or NULL")
  }
  if (!is.numeric(PI) || length(PI) != 1 || PI < 0 || PI > 1) {
    stop("PI must be a single numeric value in [0, 1]")
  }
  if (!is.numeric(N_vec) || any(N_vec <= 0) || any(N_vec != as.integer(N_vec))) {
    stop("N_vec must be a vector of positive integers")
  }
  if (!is.numeric(J) || length(J) != 1 || J < 2 || J != as.integer(J)) {
    stop("J must be a single integer >= 2")
  }
  if (!is.numeric(base_size) || length(base_size) != 1 || base_size <= 0) {
    stop("base_size must be a single positive number")
  }
  if (any(f1_seq <= 0) || any(f1_seq >= 1)) {
    stop("f1_seq must contain values strictly between 0 and 1")
  }

  # ========== Data Generation ==========
  records <- vector("list", length(N_vec) * length(f1_seq) * 4)
  idx <- 1L

  for (N in N_vec) {
    for (f1 in f1_seq) {
      N1     <- round(N * f1)
      N_rest <- N - N1

      if (N1 < 1 || N_rest < (J - 1)) next

      N_other <- floor(N_rest / (J - 1))
      N_last  <- N_rest - N_other * (J - 2)
      Nj      <- c(N1, rep(N_other, J - 2), N_last)

      if (any(Nj < 1)) next

      # Formula approach
      res_f <- rcp1armMilestoneSurvival(
        lambda = lambda, t_eval = t_eval, S0 = S0, Nj = Nj,
        t_a = t_a, t_f = t_f, lambda_dropout = lambda_dropout,
        PI = PI, approach = "formula"
      )
      # Simulation approach
      res_s <- rcp1armMilestoneSurvival(
        lambda = lambda, t_eval = t_eval, S0 = S0, Nj = Nj,
        t_a = t_a, t_f = t_f, lambda_dropout = lambda_dropout,
        PI = PI, approach = "simulation", nsim = nsim, seed = seed
      )

      for (method in c("Method 1", "Method 2")) {
        rcp_f <- if (method == "Method 1") res_f$Method1 else res_f$Method2
        rcp_s <- if (method == "Method 1") res_s$Method1 else res_s$Method2

        records[[idx]]     <- data.frame(f1 = f1, N = N, Method = method,
                                         Approach = "Formula",    RCP = rcp_f,
                                         stringsAsFactors = FALSE)
        records[[idx + 1]] <- data.frame(f1 = f1, N = N, Method = method,
                                         Approach = "Simulation", RCP = rcp_s,
                                         stringsAsFactors = FALSE)
        idx <- idx + 2L
      }
    }
  }

  df <- do.call(rbind, records[seq_len(idx - 1L)])
  df$Method   <- factor(df$Method,   levels = c("Method 1", "Method 2"))
  df$Approach <- factor(df$Approach, levels = c("Formula", "Simulation"))

  # Facet label: italic(N) = <value>  (parsed by label_parsed)
  df$N_label <- factor(paste0("italic(N) == ", df$N),
                       levels = paste0("italic(N) == ", sort(unique(df$N))))

  # ========== Plot ==========
  plt <- ggplot2::ggplot(df, ggplot2::aes(x = f1, y = RCP,
                                          color    = Method,
                                          linetype = Approach)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_hline(yintercept = 0.8, color = "#939597",
                        linetype = "longdash", linewidth = 1.0) +
    ggplot2::scale_color_manual(
      values = c("Method 1" = "#004C97", "Method 2" = "#F0B323")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Formula" = "solid", "Simulation" = "dashed")
    ) +
    ggplot2::facet_wrap(
      ~ N_label,
      nrow     = 1,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(
      title = substitute(
        paste("Milestone Survival Endpoint  (",
              lambda == a, ", ", italic(t)[eval] == b, ", ",
              italic(S)[0] == c, ", ", italic(t)[a] == d, ", ",
              italic(t)[f] == e, ", ", pi == g, ", ",
              italic(J) == h, ")"),
        list(a = round(lambda, 4), b = t_eval, c = round(S0, 4),
             d = t_a, e = t_f, g = PI, h = J)
      ),
      x = expression(italic(f)[1]),
      y = "Regional consistency probability"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.text.x     = ggplot2::element_text(size = base_size),
      strip.text.y     = ggplot2::element_text(size = base_size),
      text             = ggplot2::element_text(size = base_size),
      panel.spacing    = ggplot2::unit(-0.1, "lines"),
      legend.key.width = ggplot2::unit(2, "cm"),
      legend.text      = ggplot2::element_text(size = base_size),
      legend.title     = ggplot2::element_blank(),
      legend.position  = "bottom",
      legend.box       = "vertical"
    ) +
    ggplot2::guides(
      color    = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 2)
    )

  return(plt)
}
