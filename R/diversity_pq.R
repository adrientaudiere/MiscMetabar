################################################################################
#' Compute alpha- or beta-diversity metrics on a phyloseq object via ecodive
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Wrapper around [ecodive::beta_div()] and [ecodive::alpha_div()] providing
#' fast, C-accelerated alpha- and beta-diversity calculations directly on a
#' `phyloseq` object, plus support for `"robust_aitchison"` (alias
#' `"robust.aitchison"`).
#'
#' @details
#' **Robust Aitchison distance (`metric = "robust_aitchison"`):**
#' Two definitions exist in the microbiome literature (both citing Martino et
#' al. 2019):
#' - **`impute = FALSE` (default, `ecodive` definition):** applies the robust
#'   centered log-ratio (rCLR) transformation to strictly positive counts in
#'   each sample (`log(x) - mean(log(x[x > 0]))` for `x > 0`, and `0`
#'   otherwise) and computes the Euclidean distance via `ecodive`. This avoids
#'   pseudocounts and runs in milliseconds in C, without iterative matrix
#'   completion.
#' - **`impute = TRUE` (DEICODE / `vegan` OptSpace definition):** applies rCLR
#'   leaving zeros as `NA`, completes the low-rank matrix with OptSpace
#'   ([vegan::decostand()] with `method = "rclr", impute = TRUE`, called by
#'   `vegan::vegdist(..., method = "robust.aitchison")`), and computes the
#'   Euclidean distance on the completed matrix. Because `vegan` already
#'   implements OptSpace matrix completion natively, no external `rpca`
#'   dependency is needed.
#'
#' @inheritParams clean_pq
#' @param metric (character, default `"bray"`) Name(s) of the diversity
#'   metric(s) to compute. See `ecodive::list_metrics()` for all supported
#'   alpha and beta metric IDs, plus `"robust_aitchison"` (or
#'   `"robust.aitchison"`). All elements of `metric` must belong to the same
#'   family (either all `"beta"` or all `"alpha"`).
#' @param norm (character, default `"none"`) Normalization passed to
#'   [ecodive::beta_div()] (`"none"`, `"percent"`, `"binary"`, `"clr"`). For
#'   alpha-diversity metrics that require proportions, `"none"` is mapped to
#'   `ecodive::alpha_div()`'s default `"percent"`.
#' @param pseudocount (numeric or `NULL`, default `NULL`) Pseudocount added to
#'   counts to handle zeros when `metric = "aitchison"` or `norm = "clr"`.
#' @param impute (logical, default `FALSE`) Only used when `metric` is
#'   `"robust_aitchison"` or `"robust.aitchison"`. If `FALSE` (default), uses
#'   `ecodive`'s zero-filled rCLR Euclidean distance. If `TRUE`, uses
#'   `vegan::vegdist(..., method = "robust.aitchison")` with OptSpace low-rank
#'   matrix completion.
#' @param add_to_phyloseq (logical, default `FALSE`) Only used for
#'   alpha-diversity metrics. If `TRUE`, appends the computed alpha-diversity
#'   columns to `physeq@sam_data` and returns the updated `phyloseq` object.
#' @param nproc (integer, default `1L`) Number of CPU threads passed to
#'   `ecodive` (`cpus`). Defaults to `1L` to comply with CRAN check limits.
#' @param ... Additional arguments forwarded to [ecodive::beta_div()],
#'   [ecodive::alpha_div()], or [vegan::vegdist()] (when `impute = TRUE`).
#'
#' @return
#' - For **beta-diversity**: a [stats::dist] object when `length(metric) == 1`,
#'   or a named list of [stats::dist] objects when `length(metric) > 1`.
#' - For **alpha-diversity**: a [tibble::tibble] with column `Sample` and one
#'   column per metric (when `add_to_phyloseq = FALSE`), or a `phyloseq`
#'   object with those columns added to `sample_data(physeq)` (when
#'   `add_to_phyloseq = TRUE`).
#'
#' @export
#' @author Adrien Taudière
#' @references
#' Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A, Knight R,
#' Zengler K (2019). A novel sparse compositional technique reveals microbial
#' perturbations. *mSystems* 4(1): e00016-19.
#' \doi{10.1128/mSystems.00016-19}
#' @seealso [ecodive::beta_div()], [ecodive::alpha_div()],
#'   [vegan::vegdist()], [hill_pq()]
#' @examples
#' \donttest{
#' if (requireNamespace("ecodive", quietly = TRUE)) {
#'   pq_5 <- prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)
#'   # Beta-diversity (Bray-Curtis and Robust Aitchison)
#'   diversity_pq(pq_5, metric = "bray")
#'   diversity_pq(pq_5, metric = "robust_aitchison")
#'
#'   # Alpha-diversity
#'   diversity_pq(pq_5, metric = c("shannon", "simpson", "observed"))
#' }
#' }
diversity_pq <- function(
  physeq,
  metric = "bray",
  norm = "none",
  pseudocount = NULL,
  impute = FALSE,
  add_to_phyloseq = FALSE,
  nproc = 1L,
  ...
) {
  if (!requireNamespace("ecodive", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg ecodive} is required for {.fn diversity_pq}. Install it with {.code install.packages('ecodive')}."
    )
  }
  verify_pq(physeq)
  if (!is.character(metric) || length(metric) < 1 || anyNA(metric)) {
    cli::cli_abort("{.arg metric} must be a non-empty character vector.")
  }

  metrics_df <- ecodive::list_metrics()
  alpha_ids <- metrics_df$id[metrics_df$div == "alpha"]
  beta_ids <- unique(c(
    metrics_df$id[metrics_df$div == "beta"],
    "robust_aitchison",
    "robust.aitchison"
  ))

  norm_metric <- ifelse(
    metric == "robust.aitchison",
    "robust_aitchison",
    metric
  )

  unknown <- setdiff(norm_metric, c(alpha_ids, beta_ids))
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown diversity {.arg metric}: {.val {unknown}}.",
      "i" = "See {.code ecodive::list_metrics()} or use {.val robust_aitchison}."
    ))
  }

  is_alpha <- norm_metric %in% alpha_ids
  is_beta <- norm_metric %in% beta_ids
  if (any(is_alpha) && any(is_beta)) {
    cli::cli_abort(
      "Cannot mix alpha-diversity ({.val {metric[is_alpha]}}) and beta-diversity ({.val {metric[is_beta]}}) metrics in a single {.fn diversity_pq} call."
    )
  }

  cpus <- as.integer(nproc)

  if (all(is_beta)) {
    if (add_to_phyloseq) {
      cli::cli_abort(
        "{.arg add_to_phyloseq} = TRUE is only supported for alpha-diversity metrics."
      )
    }
    phylo_metrics <- metrics_df$id[
      metrics_df$div == "beta" & metrics_df$phylo
    ]
    if (any(norm_metric %in% phylo_metrics) && is.null(physeq@phy_tree)) {
      cli::cli_abort(
        "Phylogenetic beta-diversity metric {.val {metric[norm_metric %in% phylo_metrics]}} requires a {.field phy_tree} slot in {.arg physeq}."
      )
    }

    compute_one_beta <- function(m) {
      if (m == "robust_aitchison") {
        if (impute) {
          if (!requireNamespace("vegan", quietly = TRUE)) {
            cli::cli_abort(
              "Package {.pkg vegan} is required when {.code impute = TRUE}."
            )
          }
          otu <- as(otu_table(taxa_as_columns(physeq)), "matrix")
          mode(otu) <- "numeric"
          return(vegan::vegdist(otu, method = "robust.aitchison", ...))
        }
        if ("robust_aitchison" %in% metrics_df$id) {
          return(ecodive::beta_div(
            counts = physeq,
            metric = "robust_aitchison",
            cpus = cpus,
            ...
          ))
        }
        otu <- as(otu_table(taxa_as_columns(physeq)), "matrix")
        mode(otu) <- "numeric"
        pos <- otu > 0
        log_otu <- matrix(
          0,
          nrow = nrow(otu),
          ncol = ncol(otu),
          dimnames = dimnames(otu)
        )
        log_otu[pos] <- log(otu[pos])
        row_means <- rowSums(log_otu) / pmax(rowSums(pos), 1L)
        rclr_otu <- log_otu - row_means * pos
        d <- ecodive::euclidean(rclr_otu, margin = 1L, cpus = cpus, ...)
        attr(d, "method") <- "robust_aitchison"
        return(d)
      }

      ecodive::beta_div(
        counts = physeq,
        metric = m,
        norm = norm,
        pseudocount = pseudocount,
        cpus = cpus,
        ...
      )
    }

    if (length(norm_metric) == 1L) {
      return(compute_one_beta(norm_metric))
    }
    res_list <- lapply(norm_metric, compute_one_beta)
    names(res_list) <- metric
    return(res_list)
  }

  # Alpha-diversity path
  phylo_alpha <- metrics_df$id[
    metrics_df$div == "alpha" & metrics_df$phylo
  ]
  if (any(norm_metric %in% phylo_alpha) && is.null(physeq@phy_tree)) {
    cli::cli_abort(
      "Phylogenetic alpha-diversity metric {.val {metric[norm_metric %in% phylo_alpha]}} requires a {.field phy_tree} slot in {.arg physeq}."
    )
  }
  alpha_norm <- if (identical(norm, "none")) "percent" else norm

  cols <- lapply(norm_metric, function(m) {
    as.numeric(ecodive::alpha_div(
      counts = physeq,
      metric = m,
      norm = alpha_norm,
      cpus = cpus,
      ...
    ))
  })
  names(cols) <- metric
  res_tbl <- tibble::as_tibble(c(
    list(Sample = sample_names(physeq)),
    cols
  ))

  if (add_to_phyloseq) {
    sam_df <- as.data.frame(physeq@sam_data)
    for (nm in metric) {
      sam_df[[nm]] <- res_tbl[[nm]]
    }
    physeq@sam_data <- sample_data(sam_df)
    return(physeq)
  }

  res_tbl
}
