#' Plot Regional Consistency Probability for Single-Arm MRCT (Count Endpoint)
#'
#' @description
#' Generate a faceted plot of Regional Consistency Probability (RCP) as a function
#' of the regional allocation proportion \eqn{f_1} for count endpoints
#' (negative binomial model).
#' Formula and simulation results are shown for Method 1 (log-RR and linear-RR
#' scales) and Method 2. Facet rows correspond to the two Method 1 scales
#' (\eqn{\log(RR)} and \eqn{1 - RR}), and facet columns correspond to total
#' sample sizes specified in \code{N_vec}.
#'
#' Regional sample sizes are allocated as:
#' \eqn{N_{j1} = \lfloor N \times f_1 \rfloor} and
#' \eqn{N_{j2} = \cdots = N_{jJ} = (N - N_{j1}) / (J - 1)}.
#'
#' @param lambda Numeric scalar. Expected count per patient under the alternative
#'   hypothesis. Must be positive. Default is \code{2}.
#' @param lambda0 Numeric scalar. Expected count per patient under the historical
#'   control. Must be positive. Default is \code{3}.
#' @param dispersion Numeric scalar. Dispersion parameter of the negative binomial
#'   distribution. Must be positive. Default is \code{1}.
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
#'   scale_linetype_manual scale_x_continuous scale_y_continuous facet_grid
#'   labs theme_bw theme element_text element_blank unit guide_legend guides
#'   label_parsed as_labeller vars
#'
#' @examples
#' p <- plot_rcp1armCount(
#'   lambda     = 2,
#'   lambda0    = 3,
#'   dispersion = 1,
#'   PI         = 0.5,
#'   N_vec      = c(20, 40, 100),
#'   J          = 3
#' )
#' print(p)
#'
#' @export
plot_rcp1armCount <- function(lambda     = 2,
                              lambda0    = 3,
                              dispersion = 1,
                              PI         = 0.5,
                              N_vec      = c(20, 40, 100),
                              J          = 3,
                              f1_seq     = seq(0.1, 0.9, by = 0.1),
                              nsim       = 1e4,
                              seed       = 1,
                              base_size  = 28) {

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
  # Scale labels for facet rows
  scale_levels <- c("log(RR)", "1 - RR")
  records <- vector("list", length(N_vec) * length(f1_seq) * 6)
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
      res_f <- rcp1armCount(
        lambda = lambda, lambda0 = lambda0, dispersion = dispersion,
        Nj = Nj, PI = PI, approach = "formula"
      )
      # Simulation approach
      res_s <- rcp1armCount(
        lambda = lambda, lambda0 = lambda0, dispersion = dispersion,
        Nj = Nj, PI = PI, approach = "simulation", nsim = nsim, seed = seed
      )

      # Method 1 (log-RR), Method 1 (1-RR), Method 2
      rcp_vals <- list(
        list(method = "Method 1", scale = "log(RR)",
             rcp_f = res_f$Method1_logRR,    rcp_s = res_s$Method1_logRR),
        list(method = "Method 1", scale = "1 - RR",
             rcp_f = res_f$Method1_linearRR, rcp_s = res_s$Method1_linearRR),
        list(method = "Method 2", scale = "log(RR)",
             rcp_f = res_f$Method2,          rcp_s = res_s$Method2),
        list(method = "Method 2", scale = "1 - RR",
             rcp_f = res_f$Method2,          rcp_s = res_s$Method2)
      )

      for (rv in rcp_vals) {
        records[[idx]]     <- data.frame(f1 = f1, N = N,
                                         Method = rv$method, Scale = rv$scale,
                                         Approach = "Formula",    RCP = rv$rcp_f,
                                         stringsAsFactors = FALSE)
        records[[idx + 1]] <- data.frame(f1 = f1, N = N,
                                         Method = rv$method, Scale = rv$scale,
                                         Approach = "Simulation", RCP = rv$rcp_s,
                                         stringsAsFactors = FALSE)
        idx <- idx + 2L
      }
    }
  }

  df <- do.call(rbind, records[seq_len(idx - 1L)])
  df$Method   <- factor(df$Method,   levels = c("Method 1", "Method 2"))
  df$Approach <- factor(df$Approach, levels = c("Formula", "Simulation"))
  df$Scale    <- factor(df$Scale,    levels = scale_levels)

  # Facet column label: italic(N) = <value>
  df$N_label <- factor(paste0("italic(N) == ", df$N),
                       levels = paste0("italic(N) == ", sort(unique(df$N))))

  # Facet row label (plotmath): log(RR) and 1 - RR
  scale_labeller <- ggplot2::as_labeller(
    c("log(RR)" = "log(RR)",
      "1 - RR"  = "1 - RR"),
    default = ggplot2::label_parsed
  )

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
    ggplot2::facet_grid(
      rows     = ggplot2::vars(Scale),
      cols     = ggplot2::vars(N_label),
      labeller = ggplot2::labeller(
        Scale   = scale_labeller,
        N_label = ggplot2::label_parsed
      )
    ) +
    ggplot2::scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(
      title = substitute(
        paste("Count Endpoint  (",
              lambda == a, ", ", lambda[0] == b, ", ",
              phi == c, ", ", pi == d, ", ",
              italic(J) == e, ")"),
        list(a = lambda, b = lambda0, c = dispersion, d = PI, e = J)
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
