utils::globalVariables(c("x_value", "bin", "y_stack"))

################################################################################
#' Wheat plot of a numeric distribution
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Draw a "wheat plot" (a.k.a. wheat-ear or stacked-dot plot) of a single
#' numeric variable. Values are binned along the x-axis and, within each
#' bin, individual observations are stacked vertically as points. The result
#' is a hybrid between a histogram and a dot plot that keeps every
#' observation visible, which is useful to inspect the distribution of
#' per-taxon or per-sample quantities (e.g. `taxa_sums()` or
#' `sample_sums()`).
#'
#' @param data (data.frame, required) A data frame containing the variable
#'   to plot.
#' @param xvar (required) The (unquoted) name of the numeric column of
#'   `data` to plot.
#' @param binwidth (numeric, default: NULL) Width of the bins. When NULL,
#'   the Freedman-Diaconis rule is used to choose a value automatically.
#' @param fill (character, default: "steelblue") Colour of the points.
#' @param point_size (numeric, default: 2) Size of the points.
#' @param xlab (character, default: NULL) x-axis label. Defaults to the name
#'   of `xvar`.
#' @param ylab (character, default: "Count") y-axis label.
#' @param title (character, default: "Wheat Plot") Plot title.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @author Adrien Taudière
#' @export
#'
#' @examples
#' set.seed(42)
#' wheat_plot(
#'   data.frame(value = rnorm(200, mean = 50, sd = 10)),
#'   value,
#'   binwidth = 2
#' )
#' \donttest{
#' wheat_plot(
#'   data.frame(value = taxa_sums(data_fungi_mini)),
#'   value,
#'   binwidth = 2000
#' )
#' }
wheat_plot <- function(
  data,
  xvar,
  binwidth = NULL,
  fill = "steelblue",
  point_size = 2,
  xlab = NULL,
  ylab = "Count",
  title = "Wheat Plot"
) {
  lifecycle::deprecate_soft(
    "0.17.0",
    "wheat_plot()",
    "ggplotpq::wheat_plot()"
  )
  xvar_sym <- rlang::ensym(xvar)

  dat <- data |>
    mutate(x_value = !!xvar_sym)

  if (!is.numeric(dat$x_value)) {
    cli::cli_abort("{.arg xvar} must refer to a numeric column.")
  }

  if (is.null(binwidth)) {
    iqr <- stats::IQR(dat$x_value, na.rm = TRUE)
    n <- sum(!is.na(dat$x_value))
    binwidth <- 2 * iqr / n^(1 / 3)
    if (!is.finite(binwidth) || binwidth <= 0) {
      binwidth <- diff(range(dat$x_value, na.rm = TRUE)) / 30
    }
  }

  breaks <- seq(
    floor(min(dat$x_value, na.rm = TRUE)),
    ceiling(max(dat$x_value, na.rm = TRUE)) + binwidth,
    by = binwidth
  )

  dat <- dat |>
    mutate(
      bin = cut(x_value, breaks = breaks, include.lowest = TRUE, right = FALSE)
    ) |>
    group_by(bin) |>
    arrange(x_value, .by_group = TRUE) |>
    mutate(y_stack = dplyr::row_number()) |>
    ungroup()

  ggplot(dat, aes(x = x_value, y = y_stack)) +
    geom_point(color = fill, size = point_size) +
    labs(
      x = xlab %||% rlang::as_label(xvar_sym),
      y = ylab,
      title = title
    ) +
    theme_minimal()
}
