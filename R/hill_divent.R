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
  otu_t <- t(otu)
  # divent::profile_hill expects: first column = species names,
  # subsequent columns = one community (sample) per column
  otu_df <- data.frame(
    species = rownames(otu_t),
    otu_t,
    check.names = FALSE
  )
  result <- divent::profile_hill(otu_df, orders = orders, ...)
  ggplot2::autoplot(result)
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
#' (generalised species accumulation curves) from a `phyloseq` object, and
#' returns a ggplot2 object via [ggplot2::autoplot()].
#'
#' @inheritParams clean_pq
#' @param q (numeric, default 1) Hill diversity order.
#' @param merge_sample_by (character or NULL) If not NULL, merge samples
#'   using [merge_samples2()] before computing the curve.
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
  otu_t <- t(otu)
  # divent::accum_hill expects: first column = species names,
  # subsequent columns = one community (sample) per column
  otu_df <- data.frame(
    species = rownames(otu_t),
    otu_t,
    check.names = FALSE
  )
  result <- divent::accum_hill(otu_df, q = q, ...)
  ggplot2::autoplot(result)
}
