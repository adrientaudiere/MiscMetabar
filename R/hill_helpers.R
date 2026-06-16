################################################################################
#' Compute Hill diversity numbers for all samples in an OTU table
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Iterates over all samples in an OTU table and computes Hill diversity
#' numbers using [divent::div_hill()].
#'
#' @param comm (data.frame or matrix) OTU table with samples as rows and
#'   taxa as columns.
#' @param q (numeric vector) Hill diversity orders to compute. Hill numbers are
#'   more appropriate in DNA metabarcoding studies when `q > 0` (Alberdi &
#'   Gilbert, 2019; Calderón-Sanou et al., 2019).
#' @param ... Additional arguments passed to [divent::div_hill()] (e.g.
#'   `estimator = "naive"` to reproduce vegan-equivalent results).
#'
#' @return A data.frame with one row per sample and one column per value in
#'   `q`. Column names are the string representation of the `q` values.
#'   Row names match the input row names.
#'
#' @seealso [divent::div_hill()], [hill_pq()], [hill_tuckey_pq()]
#' @references
#' Alberdi, A., & Gilbert, M. T. P. (2019). A guide to the application of
#'   Hill numbers to DNA-based diversity analyses. *Molecular Ecology Resources*.
#'   \doi{10.1111/1755-0998.13014}
#'
#' Calderón-Sanou, I., Münkemüller, T., Boyer, F., Zinger, L., & Thuiller, W.
#'   (2019). From environmental DNA sequences to ecological conclusions: How
#'   strong is the influence of methodological choices? *Journal of Biogeography*,
#'   47. \doi{10.1111/jbi.13681}
#' @export
#' @examples
#' data("data_fungi_mini", package = "MiscMetabar")
#' data_f <- prune_samples(
#'   sample_names(data_fungi_mini)[1:5],
#'   data_fungi_mini
#' )
#' otu <- as.data.frame(phyloseq::otu_table(
#'   taxa_as_columns(data_f)
#' ))
#' divent_hill_matrix_pq(otu, q = c(0, 1, 2))
divent_hill_matrix_pq <- function(comm, q, ...) {
  comm_mat <- as.matrix(comm)
  result_raw <- vapply(
    seq_len(nrow(comm_mat)),
    function(i) {
      x <- comm_mat[i, ]
      x <- x[x > 0]
      if (length(x) == 0) {
        return(rep(NA_real_, length(q)))
      }
      vapply(
        q,
        function(q_val) {
          divent::div_hill(x, q = q_val, as_numeric = TRUE, check_arguments = FALSE, ...)
        },
        numeric(1)
      )
    },
    numeric(length(q))
  )
  if (length(q) == 1L) {
    result <- matrix(result_raw, ncol = 1L)
  } else {
    result <- t(result_raw)
  }
  out <- as.data.frame(result)
  colnames(out) <- as.character(q)
  rownames(out) <- rownames(comm_mat)
  out
}

# Internal helper: compute a single diversity index for a single sample.
# Routes "shannon"/"simpson" to divent; other indices fall back to vegan.
.divent_diversity <- function(x, index, ...) {
  x <- as.numeric(x)
  x <- x[x > 0]
  switch(index,
    "shannon" = divent::ent_shannon(x, as_numeric = TRUE, check_arguments = FALSE, ...),
    "simpson" = divent::ent_simpson(x, as_numeric = TRUE, check_arguments = FALSE, ...),
    stop("Unsupported index for divent: ", index, call. = FALSE)
  )
}

# Internal helper: compute diversity indices for a full OTU table
# (samples as rows). Routes "shannon"/"simpson" to divent, others to vegan.
# Returns a named list of numeric vectors, one per index.
.compute_diversity_index <- function(comm, indices, ...) {
  divent_indices <- c("shannon", "simpson")
  lapply(stats::setNames(indices, indices), function(idx) {
    if (idx %in% divent_indices) {
      vapply(
        seq_len(nrow(comm)),
        function(i) .divent_diversity(comm[i, ], idx, ...),
        numeric(1)
      )
    } else {
      as.numeric(vegan::diversity(comm, index = idx))
    }
  })
}
