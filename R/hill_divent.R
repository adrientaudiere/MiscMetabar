################################################################################
#' Hill diversity profile for a phyloseq object
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#' alt="lifecycle-experimental"></a>
#'
#' Wraps [divent::profile_hill()] to compute a Hill diversity profile across
#' diversity orders for each sample in a `phyloseq` object, and returns a
#' ggplot2 object via [ggplot2::autoplot()].
#'
#' @inheritParams clean_pq
#' @param orders (numeric vector) Hill diversity orders to compute. Default
#'   `seq(0, 2, 0.1)`.
#' @param merge_sample_by (character or NULL) If not NULL, merge samples
#'   using [merge_samples2()] before computing profiles.
#' @param ... Additional arguments passed to [divent::profile_hill()].
#'
#' @return A ggplot2 object.
#' @seealso [divent::profile_hill()], [hill_curves_pq()]
#' @export
#' @examples
#' \donttest{
#' profile_hill_pq(data_fungi_mini)
#' }
profile_hill_pq <- function(
  physeq,
  orders = seq(0, 2, 0.1),
  merge_sample_by = NULL,
  ...
) {
  verify_pq(physeq)
  if (!is.null(merge_sample_by)) {
    physeq <- merge_samples2(physeq, merge_sample_by)
  }
  physeq <- taxa_as_columns(physeq)
  otu <- as.data.frame(otu_table(physeq))
  sample_names <- rownames(otu)
  results <- lapply(
    cli::cli_progress_along(sample_names, name = "Computing Hill profiles"),
    \(i) {
      abund <- as.numeric(otu[i, ])
      abund <- abund[abund > 0]
      res <- divent::profile_hill(abund, orders = orders, ...)
      res$site <- sample_names[i]
      res
    }
  )
  combined <- dplyr::bind_rows(results)
  class(combined) <- c("profile", class(combined))
  ggplot2::autoplot(combined)
}

################################################################################
#' Hill diversity accumulation curve for a phyloseq object (default: q = 1)
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#' alt="lifecycle-experimental"></a>
#'
#' Computes Hill diversity accumulation curves from a `phyloseq` object
#' and returns a ggplot2 object.
#'
#' Two types of curves are available:
#' - `type = "individual"` (default): individual-based (sequence-based)
#'   rarefaction/extrapolation curves via [divent::accum_hill()], with one
#'   curve per sample (or per merged group).
#' - `type = "sample"`: sample-based accumulation curve. Samples are pooled
#'   incrementally (over random permutations) and Hill diversity is computed
#'   at each step using [divent::div_hill()]. The x-axis represents the
#'   number of samples.
#'
#' @inheritParams clean_pq
#' @param q (numeric, default 1) Hill diversity order. Default is 1 (exponential
#'   of Shannon entropy), recommended for its robustness against rare and
#'   potentially erroneous sequences (Alberdi & Gilbert, 2019;
#'   Calderón-Sanou et al., 2019).
#' @param type (character) Type of accumulation curve. Either `"individual"`
#'   (sequence-based, one curve per sample) or `"sample"` (sample-based,
#'   one curve for the whole dataset or per group).
#' @param merge_sample_by (character or NULL) Variable name in `sam_data` used
#'   to group samples. Behaviour differs by `type`:
#'   - `type = "individual"`: samples are merged before computing curves
#'     (one curve per merged group, using [merge_samples2()]).
#'   - `type = "sample"`: samples are split by group and one accumulation
#'     curve is drawn per group, all on the same plot.
#' @param n_permutations (integer, default 100) Number of random sample
#'   orderings used to compute the mean and confidence envelope for
#'   sample-based accumulation. Ignored when `type = "individual"`.
#' @param conf_level (numeric, default 0.95) Confidence level for the
#'   envelope around sample-based accumulation curves. Ignored when
#'   `type = "individual"`.
#' @param ... Additional arguments passed to [divent::accum_hill()] (when
#'   `type = "individual"`) or [divent::div_hill()] (when `type = "sample"`).
#'
#' @return A ggplot2 object.
#' @references
#' Alberdi, A., & Gilbert, M. T. P. (2019). A guide to the application of
#'   Hill numbers to DNA-based diversity analyses. *Molecular Ecology Resources*.
#'   \doi{10.1111/1755-0998.13014}
#'
#' Calderón-Sanou, I., Münkemüller, T., Boyer, F., Zinger, L., & Thuiller, W.
#'   (2019). From environmental DNA sequences to ecological conclusions: How
#'   strong is the influence of methodological choices? *Journal of Biogeography*,
#'   47. \doi{10.1111/jbi.13681}
#' @seealso [divent::accum_hill()], [divent::div_hill()], [hill_curves_pq()]
#' @export
#' @examples
#' \donttest{
#'  # Individual (sequence-based) accumulation curves
#'   hill_acc_pq(rarefy_even_depth (data_fungi_mini, sample.size = 500)) + no_legend()
#'   hill_acc_pq(rarefy_even_depth(data_fungi_mini, sample.size = 500),
#'     n_permutations = 5,
#'     merge_sample_by = "Height"
#'   )
#' 
#'  # Sample-based accumulation curve
#'   hill_acc_pq(data_fungi_mini, type = "sample", n_permutations = 50)
#'   hill_acc_pq(data_fungi_mini, type = "sample", merge_sample_by = "Height")
#' }
hill_acc_pq <- function(
  physeq,
  q = 1,
  type = c("individual", "sample"),
  merge_sample_by = NULL,
  n_permutations = 100,
  conf_level = 0.95,
  ...
) {
  type <- match.arg(type)
  verify_pq(physeq)

  if (type == "individual") {
    if (!is.null(merge_sample_by)) {
      physeq <- merge_samples2(physeq, merge_sample_by)
    }
    physeq <- taxa_as_columns(physeq)
    otu <- as.data.frame(otu_table(physeq))
    return(.hill_acc_individual(otu, q = q, ...))
  }

  physeq <- taxa_as_columns(physeq)
  otu <- as.data.frame(otu_table(physeq))

  if (is.null(merge_sample_by)) {
    return(.hill_acc_sample(
      otu,
      q = q,
      n_permutations = n_permutations,
      conf_level = conf_level,
      ...
    ))
  }

  # Split samples by group and draw one curve per group
  groups <- as.character(phyloseq::sample_data(physeq)[[merge_sample_by]])
  group_levels <- unique(groups[!is.na(groups)])
  results <- lapply(group_levels, function(g) {
    idx <- which(groups == g)
    accum <- .hill_acc_sample(
      otu[idx, , drop = FALSE],
      q = q,
      n_permutations = n_permutations,
      conf_level = conf_level,
      ...,
      .return_df = TRUE
    )
    accum$group <- g
    accum
  })
  combined <- dplyr::bind_rows(results)

  ggplot2::ggplot(
    combined,
    ggplot2::aes(
      x = .data$level,
      y = .data$diversity,
      color = .data$group,
      fill = .data$group
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lwr, ymax = .data$upr),
      alpha = 0.2,
      color = NA
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = paste("Sample-based Hill accumulation curve (q =", q, ")"),
      x = "Number of samples",
      y = paste("Hill diversity (q =", q, ")"),
      color = merge_sample_by,
      fill = merge_sample_by
    )
}

#' @noRd
.hill_acc_individual <- function(otu, q, ...) {
  sample_names <- rownames(otu)
  results <- lapply(
    cli::cli_progress_along(
      sample_names,
      name = "Computing individual-based accumulation curves"
    ),
    \(i) {
      abund <- as.numeric(otu[i, ])
      abund <- abund[abund > 0]
      res <- divent::accum_hill(abund, q = q, ...)
      res$site <- sample_names[i]
      res
    }
  )
  combined <- dplyr::bind_rows(results)
  class(combined) <- c("accumulation", class(combined))
  ggplot2::autoplot(combined) +
    ggplot2::labs(
      title = paste("Hill diversity accumulation curves (q =", q, ")"),
      x = "Number of sequences",
      y = paste("Hill diversity (q =", q, ")")
    )
}

#' @noRd
.hill_acc_sample <- function(
  otu,
  q,
  n_permutations,
  conf_level,
  ...,
  .return_df = FALSE
) {
  n_samples <- nrow(otu)
  otu_mat <- as.matrix(otu)
  alpha <- 1 - conf_level

  div_matrix <- matrix(NA_real_, nrow = n_permutations, ncol = n_samples)
  for (p in cli::cli_progress_along(
    seq_len(n_permutations),
    name = "Computing sample-based accumulation curves"
  )) {
    perm <- sample.int(n_samples)
    for (k in seq_len(n_samples)) {
      pooled <- colSums(otu_mat[perm[seq_len(k)], , drop = FALSE])
      pooled <- pooled[pooled > 0]
      div_matrix[p, k] <- divent::div_hill(
        pooled,
        q = q,
        as_numeric = TRUE,
        estimator = "naive",
        ...
      )
    }
  }

  accum_df <- data.frame(
    level = seq_len(n_samples),
    diversity = colMeans(div_matrix),
    lwr = apply(div_matrix, 2, stats::quantile, probs = alpha / 2),
    upr = apply(div_matrix, 2, stats::quantile, probs = 1 - alpha / 2)
  )

  if (.return_df) {
    return(accum_df)
  }

  ggplot2::ggplot(
    accum_df,
    ggplot2::aes(x = .data$level, y = .data$diversity)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lwr, ymax = .data$upr),
      alpha = 0.2
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = paste("Sample-based Hill accumulation curve (q =", q, ")"),
      x = "Number of samples",
      y = paste("Hill diversity (q =", q, ")")
    )
}
