################################################################################
#' Unified dispatcher for all OTU-table transformations and normalisations
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Single entry-point for all count-table transformations available in
#'   MiscMetabar.  Ecological methods (`"tss"`, `"hellinger"`, `"clr"`,
#'   `"rclr"`, `"log1p"`, `"z"`, `"pa"`, `"rank"`) are delegated to
#'   [vegan::decostand()].  Library-size normalisation methods (`"rarefy"`,
#'   `"srs"`, `"gmpr"`, `"css"`, `"tmm"`, `"vst"`) and the McKnight
#'   log-log residual method (`"mcknight_residuals"`) are delegated to
#'   their dedicated `*_pq()` functions.  All `...` arguments are forwarded
#'   to the underlying function.
#'
#' @inheritParams clean_pq
#' @param method (character, default `"tss"`) transformation to apply. One of:
#'   \describe{
#'     \item{`"tss"`}{Total Sum Scaling — divide by library size.}
#'     \item{`"hellinger"`}{Square-root of proportions. Good for ordination.}
#'     \item{`"clr"`}{Centred log-ratio (adds `pseudocount` to handle zeros).}
#'     \item{`"rclr"`}{Robust CLR (adds `pseudocount` to handle zeros).}
#'     \item{`"log1p"`}{\eqn{\log(1 + x)} transformation.}
#'     \item{`"z"`}{Per-taxon z-score standardisation.}
#'     \item{`"pa"`}{Presence/absence (0/1).}
#'     \item{`"rank"`}{Replace counts by within-sample ranks.}
#'     \item{`"normalize_prop"`}{TSS × constant + log, via [normalize_prop_pq()].}
#'     \item{`"rarefy"`}{Rarefaction to equal depth, via [rarefy_pq()].}
#'     \item{`"srs"`}{Scaling with Ranked Subsampling, via [srs_pq()].
#'       Requires the \pkg{SRS} package.}
#'     \item{`"gmpr"`}{Geometric Mean of Pairwise Ratios, via [gmpr_pq()].}
#'     \item{`"css"`}{Cumulative Sum Scaling, via [css_pq()].
#'       Requires the \pkg{metagenomeSeq} package.}
#'     \item{`"tmm"`}{Trimmed Mean of M-values, via [tmm_pq()].
#'       Requires the \pkg{edgeR} package.}
#'     \item{`"vst"`}{Variance Stabilising Transformation, via [vst_pq()].
#'       Requires the \pkg{DESeq2} package.}
#'     \item{`"mcknight_residuals"`}{Log-log depth residuals added to
#'       `sample_data`, via [mcknight_residuals_pq()].}
#'   }
#' @param pseudocount (numeric, default `1`) added before `"clr"` / `"rclr"`
#'   to avoid non-positive values.  Ignored for all other methods.
#' @param ... Additional arguments forwarded to the underlying function
#'   ([vegan::decostand()], [rarefy_pq()], [srs_pq()], etc.).
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with a
#'   transformed `otu_table` (and augmented `sample_data` for
#'   `"mcknight_residuals"`).
#' @export
#' @author Adrien Taudière
#' @seealso [normalize_prop_pq()], [rarefy_pq()], [srs_pq()], [gmpr_pq()],
#'   [css_pq()], [tmm_pq()], [vst_pq()], [mcknight_residuals_pq()],
#'   [as_binary_otu_table()], [vegan::decostand()]
#' @examplesIf rlang::is_installed("vegan")
#' data_f_tss  <- transform_pq(data_fungi_mini, method = "tss")
#' sample_sums(data_f_tss)
#'
#' data_f_hell <- transform_pq(data_fungi_mini, method = "hellinger")
#' data_f_clr  <- transform_pq(data_fungi_mini, method = "clr")
#' data_f_pa   <- transform_pq(data_fungi_mini, method = "pa")
#' data_f_rar  <- transform_pq(data_fungi_mini, method = "rarefy", seed = 1)
#' data_f_gmpr <- transform_pq(data_fungi_mini, method = "gmpr")
transform_pq <- function(
  physeq,
  method = c(
    "tss",
    "hellinger",
    "clr",
    "rclr",
    "log1p",
    "z",
    "pa",
    "rank",
    "normalize_prop",
    "rarefy",
    "srs",
    "gmpr",
    "css",
    "tmm",
    "vst",
    "mcknight_residuals"
  ),
  pseudocount = 1,
  ...
) {
  verify_pq(physeq)
  method <- match.arg(method)

  # ── delegated methods ────────────────────────────────────────────────────
  if (method == "normalize_prop") {
    return(normalize_prop_pq(physeq, ...))
  }
  if (method == "rarefy") {
    return(rarefy_pq(physeq, ...))
  }
  if (method == "srs") {
    return(srs_pq(physeq, ...))
  }
  if (method == "gmpr") {
    return(gmpr_pq(physeq, ...))
  }
  if (method == "css") {
    return(css_pq(physeq, ...))
  }
  if (method == "tmm") {
    return(tmm_pq(physeq, ...))
  }
  if (method == "vst") {
    return(vst_pq(physeq, ...))
  }
  if (method == "mcknight_residuals") {
    return(mcknight_residuals_pq(physeq, ...))
  }

  # ── vegan::decostand methods ─────────────────────────────────────────────
  otu <- as(physeq@otu_table, "matrix")
  if (taxa_are_rows(physeq)) {
    otu <- t(otu)
  }

  new_otu <- switch(
    method,
    tss = vegan::decostand(otu, method = "total", ...),
    hellinger = vegan::decostand(otu, method = "hellinger", ...),
    clr = vegan::decostand(otu + pseudocount, method = "clr", ...),
    rclr = vegan::decostand(otu + pseudocount, method = "rclr", ...),
    z = vegan::decostand(otu, method = "standardize", ...),
    rank = vegan::decostand(otu, method = "rank", ...),
    pa = {
      res <- otu
      res[res > 0] <- 1
      res
    },
    log1p = log1p(otu)
  )

  new_physeq <- physeq
  if (taxa_are_rows(physeq)) {
    new_physeq@otu_table <- otu_table(t(new_otu), taxa_are_rows = TRUE)
  } else {
    new_physeq@otu_table <- otu_table(new_otu, taxa_are_rows = FALSE)
  }

  return(new_physeq)
}
################################################################################

################################################################################
#' Normalize OTU table using samples depth
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   This function implement the method proposed by
#'   McKnight et al. 2018 (\doi{doi:10.5061/dryad.tn8qs35})
#'
#' @inheritParams clean_pq
#'
#' @param base_log (integer, default 2) the base for log-transformation. If
#'   set to NULL or NA, no log-transformation is compute after normalization.
#' @param constante a constante to multiply the otu_table values
#' @param digits (default = 4) integer indicating the number of decimal places
#'   to be used (see `?round` for more information)
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with otu_table count
#'   normalize and log transformed (if base_log is an integer)
#' @export
#' @author Adrien Taudière
#' @examples
#' taxa_sums(data_fungi_mini)
#' data_f_norm <- normalize_prop_pq(data_fungi_mini)
#' taxa_sums(data_f_norm)
#' sample_sums(data_f_norm)
#' ggplot(data.frame(
#'   "norm" = scale(taxa_sums(data_f_norm)),
#'   "raw" = scale(taxa_sums(data_fungi_mini)),
#'   "name_otu" = taxa_names(data_f_norm)
#' )) +
#'   geom_point(aes(x = raw, y = norm))
#'
#' data_f_norm <- normalize_prop_pq(taxa_as_columns(data_fungi_mini))
#'
#' data_f_norm2 <- normalize_prop_pq(data_fungi_mini, base_log = NULL)
#' taxa_sums(data_f_norm2)
#' sample_sums(data_f_norm2)
normalize_prop_pq <- function(
  physeq,
  base_log = 2,
  constante = 10000,
  digits = 4
) {
  verify_pq(physeq)
  if (taxa_are_rows(physeq)) {
    new_otutab <- round(
      (apply(physeq@otu_table, 2, function(x) {
        x / sum(x)
      })) *
        constante,
      digits = digits
    )
  } else {
    new_otutab <- round(
      (apply(physeq@otu_table, 1, function(x) {
        x / sum(x)
      })) *
        constante,
      digits = digits
    )
  }

  if (!is.null(base_log) && !is.na(base_log)) {
    new_otutab <- round(log(new_otutab + 1, base = base_log), digits = digits)
  }

  if (taxa_are_rows(physeq)) {
    new_physeq <- physeq
    new_physeq@otu_table <- otu_table(new_otutab, taxa_are_rows = TRUE)
  } else {
    new_physeq <- physeq
    new_physeq@otu_table <- otu_table(t(new_otutab), taxa_are_rows = FALSE)
  }

  return(new_physeq)
}
################################################################################

################################################################################
#' Rarefy a phyloseq object, optionally averaging over repetitions
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Wrapper around [phyloseq::rarefy_even_depth()] that can repeat the
#'   rarefaction `n` times and return the averaged OTU table, reducing the
#'   stochasticity of a single subsampling pass. With `n = 1` (default) the
#'   behaviour matches a standard single rarefaction.
#'
#' @inheritParams clean_pq
#' @param sample_size (integer) the depth to rarefy to. If `NULL` (default),
#'   the minimum `sample_sums(physeq)` is used.
#' @param n (integer, default 1) number of rarefaction repetitions to average
#'   over. Values > 1 return a non-integer (averaged) OTU table.
#' @param seed (integer, default 123) random seed passed to
#'   [phyloseq::rarefy_even_depth()].
#' @param ... Additional arguments passed on to
#'   [phyloseq::rarefy_even_depth()].
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with a
#'   rarefied (or averaged-rarefied) `otu_table`.
#' @export
#' @author Adrien Taudière
#' @seealso [phyloseq::rarefy_even_depth()], [transform_pq()]
#' @examples
#' data_f_rar <- rarefy_pq(data_fungi_mini, seed = 1)
#' sample_sums(data_f_rar)
#'
#' data_f_rar5 <- rarefy_pq(data_fungi_mini, n = 5, seed = 1)
#' sample_sums(data_f_rar5)
rarefy_pq <- function(
  physeq,
  sample_size = NULL,
  n = 1,
  seed = 123,
  ...
) {
  verify_pq(physeq)
  if (is.null(sample_size)) {
    sample_size <- min(sample_sums(physeq))
  }

  rar_once <- function(s) {
    suppressMessages(phyloseq::rarefy_even_depth(
      physeq,
      sample.size = sample_size,
      rngseed = s,
      replace = FALSE,
      trimOTUs = FALSE,
      verbose = FALSE,
      ...
    ))
  }

  # Ensure .Random.seed exists before phyloseq tries to save it
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    set.seed(NULL)
  }

  if (n <= 1) {
    return(rar_once(seed))
  }

  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n)
  reps <- lapply(seeds, rar_once)

  # Average OTU tables across repetitions
  mats <- lapply(reps, function(p) as(otu_table(p), "matrix"))
  avg <- Reduce("+", mats) / n

  new_physeq <- reps[[1]]
  new_physeq@otu_table <- otu_table(
    avg,
    taxa_are_rows = taxa_are_rows(reps[[1]])
  )
  return(new_physeq)
}
################################################################################

################################################################################
# Internal helper: run a function `fun` on a taxa-in-rows matrix of counts
# and rebuild a phyloseq object with the result, preserving orientation.
.apply_otu_transform <- function(physeq, fun) {
  otu <- as(physeq@otu_table, "matrix")
  if (!taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  new_otu <- fun(otu)
  new_physeq <- physeq
  if (taxa_are_rows(physeq)) {
    new_physeq@otu_table <- otu_table(new_otu, taxa_are_rows = TRUE)
  } else {
    new_physeq@otu_table <- otu_table(t(new_otu), taxa_are_rows = FALSE)
  }
  new_physeq
}

################################################################################
#' Scaling with Ranked Subsampling (SRS) normalization of a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Wrapper around [SRS::SRS()] (Heidrich et al. 2021,
#'   \doi{10.7717/peerj.9593}) which scales all samples to a common count
#'   `Cmin` while preserving the rank order of OTU abundances.
#'
#' @inheritParams clean_pq
#' @param Cmin (integer) the common scaling depth. Defaults to
#'   `min(sample_sums(physeq))`.
#' @param ... Additional arguments passed on to [SRS::SRS()].
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with the SRS
#'   normalised `otu_table`.
#' @export
#' @author Adrien Taudière
#' @seealso [SRS::SRS()], [rarefy_pq()]
#' @examplesIf rlang::is_installed("SRS")
#' data_f_srs <- srs_pq(data_fungi_mini)
#' sample_sums(data_f_srs)
srs_pq <- function(physeq, Cmin = NULL, ...) {
  verify_pq(physeq)
  if (is.null(Cmin)) {
    Cmin <- min(sample_sums(physeq))
  }
  .apply_otu_transform(physeq, function(mat) {
    df <- as.data.frame(mat)
    res <- SRS::SRS(df, Cmin = Cmin, ...)
    rownames(res) <- rownames(mat)
    as.matrix(res)
  })
}
################################################################################

################################################################################
#' Cumulative Sum Scaling (CSS) normalization of a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Wrapper around `metagenomeSeq::cumNorm()` / `metagenomeSeq::MRcounts()`
#'   implementing Cumulative Sum Scaling (Paulson et al. 2013,
#'   \doi{10.1038/nmeth.2658}).
#'
#' @inheritParams clean_pq
#' @param log (logical, default `TRUE`) whether to return
#'   \eqn{\log_2(x + 1)} transformed counts (as recommended by the
#'   metagenomeSeq authors).
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with a CSS
#'   normalised `otu_table`.
#' @export
#' @author Adrien Taudière
#' @seealso `metagenomeSeq::cumNorm()`
#' @examplesIf rlang::is_installed("metagenomeSeq")
#' data_f_css <- css_pq(data_fungi_mini)
css_pq <- function(physeq, log = TRUE) {
  verify_pq(physeq)
  .apply_otu_transform(physeq, function(mat) {
    mrexp <- metagenomeSeq::newMRexperiment(counts = mat)
    p <- metagenomeSeq::cumNormStatFast(mrexp)
    mrexp <- metagenomeSeq::cumNorm(mrexp, p = p)
    metagenomeSeq::MRcounts(mrexp, norm = TRUE, log = log)
  })
}
################################################################################

################################################################################
#' Trimmed Mean of M-values (TMM) normalization of a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Wrapper around [edgeR::calcNormFactors()] with `method = "TMM"`
#'   (Robinson & Oshlack 2010, \doi{10.1186/gb-2010-11-3-r25}). Returns
#'   counts-per-million scaled by the TMM-derived library sizes.
#'
#' @inheritParams clean_pq
#' @param log (logical, default `FALSE`) if `TRUE`, returns `log2(cpm + 1)`.
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with a TMM
#'   normalised `otu_table`.
#' @export
#' @author Adrien Taudière
#' @seealso [edgeR::calcNormFactors()], [edgeR::cpm()]
#' @examplesIf rlang::is_installed("edgeR")
#' data_f_tmm <- tmm_pq(data_fungi_mini)
tmm_pq <- function(physeq, log = FALSE) {
  verify_pq(physeq)
  .apply_otu_transform(physeq, function(mat) {
    dge <- edgeR::DGEList(counts = mat)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    out <- edgeR::cpm(dge, log = log)
    out
  })
}
################################################################################

################################################################################
#' Variance Stabilizing Transformation of a phyloseq object (DESeq2)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Wrapper around [DESeq2::varianceStabilizingTransformation()]
#'   (Love, Huber & Anders 2014, \doi{10.1186/s13059-014-0550-8}). Counts
#'   are incremented by 1 to handle zeros before VST is applied.
#'
#' @inheritParams clean_pq
#' @param blind (logical, default `TRUE`) passed to DESeq2.
#' @param fitType (character, default `"parametric"`) passed to DESeq2.
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with a VST
#'   transformed `otu_table`.
#' @export
#' @author Adrien Taudière
#' @seealso [DESeq2::varianceStabilizingTransformation()]
#' @examplesIf rlang::is_installed("DESeq2")
#' data_f_vst <- vst_pq(data_fungi_mini)
vst_pq <- function(physeq, blind = TRUE, fitType = "parametric") {
  verify_pq(physeq)
  .apply_otu_transform(physeq, function(mat) {
    DESeq2::varianceStabilizingTransformation(
      mat + 1L,
      blind = blind,
      fitType = fitType
    )
  })
}
################################################################################

################################################################################
#' Geometric Mean of Pairwise Ratios (GMPR) normalization of a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Pure-R implementation of the Geometric Mean of Pairwise Ratios
#'   normalization (Chen et al. 2018, \doi{10.7717/peerj.4600}) tailored for
#'   zero-inflated count tables such as microbial OTU tables. Returns counts
#'   divided by the per-sample GMPR size factors.
#'
#' @inheritParams clean_pq
#' @param intersect_no (integer, default 4) minimum number of shared taxa
#'   between two samples for the pairwise ratio to be computed.
#' @param ct_min (integer, default 2) minimum count for a taxon to be
#'   considered "shared" between two samples.
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with a
#'   GMPR-normalised `otu_table`. Size factors are stored as an attribute
#'   `"gmpr_size_factors"` on the otu_table.
#' @export
#' @author Adrien Taudière
#' @references Chen L. et al. (2018) GMPR: a robust normalization method for
#'   zero-inflated count data with application to microbiome sequencing data.
#'   PeerJ 6:e4600. \doi{10.7717/peerj.4600}
#' @examples
#' data_f_gmpr <- gmpr_pq(data_fungi_mini)
#' sample_sums(data_f_gmpr)
gmpr_pq <- function(physeq, intersect_no = 4, ct_min = 2) {
  verify_pq(physeq)

  # Work on taxa-in-rows (samples in columns)
  mat <- as(physeq@otu_table, "matrix")
  if (!taxa_are_rows(physeq)) {
    mat <- t(mat)
  }

  nsamp <- ncol(mat)
  sf <- rep(NA_real_, nsamp)
  for (i in seq_len(nsamp)) {
    log_ratios <- c()
    xi <- mat[, i]
    for (j in seq_len(nsamp)) {
      if (i == j) {
        next
      }
      xj <- mat[, j]
      shared <- (xi >= ct_min) & (xj >= ct_min)
      if (sum(shared) < intersect_no) {
        next
      }
      # median of pairwise log ratios -> median ratio
      log_ratios <- c(log_ratios, median(log(xi[shared] / xj[shared])))
    }
    if (length(log_ratios) == 0) {
      sf[i] <- NA_real_
    } else {
      # GMPR: geometric mean of the per-pair median ratios
      sf[i] <- exp(mean(log_ratios))
    }
  }
  names(sf) <- colnames(mat)

  if (any(is.na(sf))) {
    warning(
      "GMPR size factors could not be computed for ",
      sum(is.na(sf)),
      " sample(s); these samples are left unscaled."
    )
    sf[is.na(sf)] <- 1
  }

  new_mat <- sweep(mat, 2, sf, "/")
  attr(new_mat, "gmpr_size_factors") <- sf

  new_physeq <- physeq
  if (taxa_are_rows(physeq)) {
    new_physeq@otu_table <- otu_table(new_mat, taxa_are_rows = TRUE)
  } else {
    new_physeq@otu_table <- otu_table(t(new_mat), taxa_are_rows = FALSE)
  }
  return(new_physeq)
}
################################################################################

################################################################################
#' Depth-robust alpha diversity residuals (McKnight / Mikryukov)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Computes the residuals from a linear regression of
#'   \eqn{\log(\text{richness})} against \eqn{\log(\text{sequencing depth})}
#'   as a depth-robust alpha diversity metric (Mikryukov et al. 2023;
#'   McKnight et al. 2018, \doi{10.5061/dryad.tn8qs35}). This avoids
#'   discarding data through rarefaction.
#'
#' @inheritParams clean_pq
#' @param add_to_sam_data (logical, default `TRUE`) if `TRUE`, a column
#'   `mcknight_residuals` is added to `sample_data(physeq)` and the
#'   augmented phyloseq object is returned. If `FALSE`, the numeric
#'   residuals vector is returned.
#'
#' @return Either a phyloseq object with an augmented `sample_data` (default)
#'   or a named numeric vector of residuals.
#' @export
#' @author Adrien Taudière
#' @examples
#' data_f_res <- mcknight_residuals_pq(data_fungi_mini)
#' head(sample_data(data_f_res)$mcknight_residuals)
mcknight_residuals_pq <- function(physeq, add_to_sam_data = TRUE) {
  verify_pq(physeq)
  depth <- sample_sums(physeq)
  otu <- as(physeq@otu_table, "matrix")
  if (!taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  richness <- colSums(otu > 0)
  fit <- stats::lm(log(richness) ~ log(depth))
  res <- stats::residuals(fit)
  names(res) <- sample_names(physeq)

  if (!add_to_sam_data) {
    return(res)
  }

  new_physeq <- physeq
  sdf <- as(sample_data(new_physeq), "data.frame")
  sdf$mcknight_residuals <- res[rownames(sdf)]
  sample_data(new_physeq) <- sample_data(sdf)
  return(new_physeq)
}
################################################################################
