################################################################################
#' Verify that grouping columns exist in the `sam_data` slot of a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Check that the column name(s) supplied as a grouping factor (`fact` or its
#' synonym `modality`) or as a two-level grouping factor (`bifactor`) are
#' actually present in the `sample_data()` (`sam_data`) slot of a `phyloseq`
#' object. For `bifactor`, the column is additionally required to have
#' **exactly two** levels.
#'
#' The goal is to raise a clear, early error listing the available columns
#' instead of letting a missing column propagate into a cryptic downstream
#' failure (e.g. `NULL` being coerced to a zero-level factor). Most functions
#' of the pqverse that take a `physeq` together with a `fact`, `modality` or
#' `bifactor` argument call this helper internally.
#'
#' Each argument accepts a character vector, so several columns can be checked
#' in a single call. `NULL` arguments are skipped.
#'
#' @inheritParams clean_pq
#' @param fact (default NULL) Name(s) of the `sam_data` column(s) used as a
#'   primary grouping factor. Checked for presence only.
#' @param bifactor (default NULL) Name(s) of the `sam_data` column(s) used as a
#'   two-level grouping factor. Checked for presence **and** for having exactly
#'   two levels.
#' @param modality (default NULL) Synonym of `fact` kept for the functions of
#'   the pqverse that name their grouping argument `modality`. Checked for
#'   presence only.
#' @param call (default `rlang::caller_env()`) The calling environment, used to
#'   point the error message at the user-facing function rather than at
#'   `verify_fact_pq()` itself.
#'
#' @return Invisibly returns `physeq`. The function is called for its side
#'   effect: it throws an informative error when a check fails.
#' @export
#' @author Adrien Taudière
#' @examples
#' # Presence check (passes silently)
#' verify_fact_pq(data_fungi_mini, fact = "Height")
#'
#' \dontrun{
#' # Missing column: error lists the available sam_data columns
#' verify_fact_pq(data_fungi_mini, fact = "Heigth")
#'
#' # bifactor with more than two levels: error (Height has 3 levels)
#' verify_fact_pq(data_fungi_mini, bifactor = "Height")
#' }
#'
#' # A genuine two-level column passes the bifactor check
#' data_2h <- subset_samples_pq(
#'   data_fungi_mini,
#'   data_fungi_mini@sam_data$Height %in% c("Low", "High")
#' )
#' verify_fact_pq(data_2h, bifactor = "Height")
verify_fact_pq <- function(
  physeq,
  fact = NULL,
  bifactor = NULL,
  modality = NULL,
  call = rlang::caller_env()
) {
  if (is.null(physeq@sam_data)) {
    cli::cli_abort(
      "{.arg physeq} has no {.field sam_data} slot, so grouping columns cannot be checked.",
      call = call
    )
  }
  avail <- colnames(physeq@sam_data)

  check_present <- function(cols, arg) {
    if (is.null(cols)) {
      return(invisible())
    }
    for (col in cols) {
      if (!col %in% avail) {
        cli::cli_abort(
          c(
            "Column {.val {col}} (argument {.arg {arg}}) was not found in the {.field sam_data} slot of {.arg physeq}.",
            "i" = "Available columns: {.val {avail}}."
          ),
          call = call
        )
      }
    }
    invisible()
  }

  check_present(fact, "fact")
  check_present(modality, "modality")
  check_present(bifactor, "bifactor")

  if (!is.null(bifactor)) {
    for (col in bifactor) {
      n_lev <- nlevels(as.factor(physeq@sam_data[[col]]))
      if (n_lev != 2) {
        levs <- levels(as.factor(physeq@sam_data[[col]]))
        cli::cli_abort(
          c(
            "The {.arg bifactor} column {.val {col}} must have exactly two levels, but has {n_lev}.",
            "i" = if (n_lev > 0) {
              "Levels found: {.val {levs}}."
            } else {
              "The column is empty or contains only missing values."
            }
          ),
          call = call
        )
      }
    }
  }

  invisible(physeq)
}
################################################################################
