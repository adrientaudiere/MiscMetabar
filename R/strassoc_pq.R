################################################################################
#' Strength of species-group associations for a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' A wrapper for the [indicspecies::strassoc()] function in the case of a
#'   `physeq` object. It computes the strength of the association between each
#'   taxon and the groups defined by `fact`, using one of the association
#'   indices provided by the `indicspecies` package (e.g. `IndVal`, the phi
#'   coefficient `r`, the specificity component `A`, or the fidelity component
#'   `B`). It complements [multipatt_pq()] which tests the statistical
#'   significance of those associations.
#'
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor in `physeq@sam_data` defining the
#'   groups of samples.
#' @param func (chr, default "IndVal.g") The association index to compute. See
#'   `?indicspecies::strassoc()` for the full list. Common values are:
#'   - `"IndVal"` / `"IndVal.g"`: Dufrêne-Legendre indicator value (the `.g`
#'     variant equalizes group sizes),
#'   - `"r"` / `"r.g"`: the point-biserial / phi coefficient of association,
#'   - `"A"`: specificity (positive predictive value) component,
#'   - `"B"`: fidelity (sensitivity) component.
#' @param nboot_ci (int, default NULL) Number of bootstrap replicates used to
#'   estimate confidence intervals. When NULL (default), no bootstrap is
#'   performed and a single estimate is returned per taxon and group.
#' @param alpha_ci (float, default 0.05) Error level for the bootstrap
#'   confidence intervals (only used when `nboot_ci` is not NULL).
#' @param ... Additional arguments passed on to [indicspecies::strassoc()].
#'
#' @return
#'   - When `nboot_ci` is NULL: a `tibble` with one row per taxon, a
#'     `taxon` column and one column per group containing the association
#'     value.
#'   - When `nboot_ci` is set: the list returned by [indicspecies::strassoc()]
#'     with the `lowerCI`, `stat` and `upperCI` matrices.
#'
#' @export
#' @seealso [multipatt_pq()], [indicspecies::strassoc()]
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `indicspecies::strassoc()` if you
#'   use this function.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("indicspecies")) {
#'   strassoc_pq(
#'     subset_samples(data_fungi_mini, !is.na(Height)),
#'     fact = "Height",
#'     func = "IndVal.g"
#'   )
#' }
#' }
#' \dontrun{
#' if (requireNamespace("indicspecies")) {
#'   strassoc_pq(
#'     subset_samples(data_fungi_mini, !is.na(Height)),
#'     fact = "Height",
#'     func = "A"
#'   )
#'   strassoc_pq(
#'     subset_samples(data_fungi_mini, !is.na(Height)),
#'     fact = "Height",
#'     func = "IndVal.g",
#'     nboot_ci = 99
#'   )
#' }
#' }
strassoc_pq <- function(
  physeq,
  fact,
  func = "IndVal.g",
  nboot_ci = NULL,
  alpha_ci = 0.05,
  ...
) {
  if (!requireNamespace("indicspecies", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg indicspecies} is required for {.fn strassoc_pq}. Please install it."
    )
  }

  verify_pq(physeq)
  verify_fact_pq(physeq, fact = fact)

  if (nlevels(as.factor(physeq@sam_data[[fact]])) < 2) {
    cli::cli_abort(
      "The factor {.val {fact}} must have at least two levels for {.fn strassoc_pq}."
    )
  }

  physeq <- taxa_as_columns(physeq)

  res <- indicspecies::strassoc(
    as.matrix(physeq@otu_table),
    physeq@sam_data[[fact]],
    func = func,
    nboot.ci = nboot_ci,
    alpha.ci = alpha_ci,
    ...
  )

  if (is.null(nboot_ci)) {
    res <- tibble::as_tibble(
      as.data.frame(res),
      rownames = "taxon"
    )
  }

  return(res)
}
################################################################################
