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
#' Hill diversity accumulation curve for a phyloseq object
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#' alt="lifecycle-experimental"></a>
#'
#' Wraps [divent::accum_hill()] to compute Hill diversity accumulation curves
#' (rarefaction curves) from a `phyloseq` object, with one curve per sample
#' (or per merged group), and returns a ggplot2 object via
#' [ggplot2::autoplot()].
#'
#' @inheritParams clean_pq
#' @param q (numeric, default 1) Hill diversity order.
#' @param merge_sample_by (character or NULL) If not NULL, merge samples
#'   using [merge_samples2()] before computing the curves.
#' @param ... Additional arguments passed to [divent::accum_hill()].
#'
#' @return A ggplot2 object.
#' @seealso [divent::accum_hill()], [hill_curves_pq()]
#' @export
#' @examples
#' \donttest{
#' hill_acc_pq(data_fungi_mini, q = 1)
#' }
hill_acc_pq <- function(
  physeq,
  q = 1,
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
    cli::cli_progress_along(
      sample_names,
      name = "Computing accumulation curves"
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
  ggplot2::autoplot(combined)
}
