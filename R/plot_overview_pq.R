################################################################################
#' Overview plot of alpha and beta diversity for a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Produce, in a single call, the main relevant graphical views of alpha- and
#' beta-diversity for a phyloseq object as a function of **one** sample
#' variable (`fact`). The set of panels adapts to the *type* of `fact`:
#'
#' - **Numeric variable** — alpha-diversity is shown as a Hill-numbers scatter
#'   plot against the variable (via [ggscatt_pq()], one panel per Hill order
#'   `q`, with the correlation statistics from \pkg{ggstatsplot}); the
#'   ordination and UMAP points are colored by a continuous gradient. Venn/UpSet
#'   panels are skipped (they require discrete groups).
#' - **Factor with 2 to `venn_max` levels** — alpha-diversity bar plots with
#'   error bars and Tukey letters (via [hill_bar_pq()]), a Venn diagram of
#'   shared taxa (via [ggvenn_pq()]), an ordination and a UMAP colored by the
#'   levels.
#' - **Factor with more than `venn_max` levels** — same as above but the Venn
#'   diagram is replaced by an UpSet plot (via [upset_pq()]), which stays
#'   legible with many sets.
#'
#' This is a convenience wrapper meant for quick exploration; for publication
#' figures, call the dedicated functions ([hill_bar_pq()], [ggscatt_pq()],
#' [plot_ordination_pq()], [ggvenn_pq()], [upset_pq()], [umap_pq()]) directly to
#' fine-tune each panel.
#'
#' @inheritParams clean_pq
#' @param fact (required) Name of a sample variable present in the `sam_data`
#'   slot of `physeq`. Drives every panel. Numeric variables trigger the
#'   gradient/scatter behavior; other variables are treated as a factor and
#'   must have at least two levels.
#' @param q (vector of integer, default `c(0, 1, 2)`) The Hill numbers orders
#'   (q = 0 richness, q = 1 Shannon, q = 2 Simpson).
#' @param add_alpha (logical, default TRUE) Add the alpha-diversity panel: a
#'   Hill-number bar plot ([hill_bar_pq()]) for a factor, or a Hill-number
#'   scatter ([ggscatt_pq()]) for a numeric variable.
#' @param add_ordination (logical, default TRUE) Add the beta-diversity
#'   ordination panel.
#' @param add_venn (logical, default TRUE) Add a Venn (or UpSet) panel of
#'   shared taxa across the levels of `fact`. Ignored when `fact` is numeric.
#' @param add_umap (logical, default TRUE) Add a UMAP panel (via [umap_pq()]).
#'   This panel can be slow on datasets with many samples. It is skipped (with
#'   a message) when the \pkg{umap} package is not installed or when there are
#'   15 samples or fewer (`umap` defaults to `n_neighbors = 15`, which requires
#'   more samples than that); set to FALSE to skip it explicitly.
#' @param venn_max (integer, default 4) Maximum number of levels for which a
#'   Venn diagram is drawn. Above this threshold an UpSet plot is used instead.
#' @param ordination_method (character, default `"PCoA"`) Ordination method
#'   passed to [plot_ordination_pq()]. `"PCoA"` is recommended as it always
#'   converges; `"NMDS"` may fail on small or sparse datasets.
#' @param dist_method (character, default `"bray"`) Distance method passed to
#'   [plot_ordination_pq()] (and ultimately to [vegan::vegdist()]).
#' @param one_plot (logical, default FALSE) If TRUE, assemble the panels into a
#'   single figure with the \pkg{patchwork} package. If FALSE (default), return
#'   a named list (one entry per panel), letting the user arrange them freely.
#'   Mirrors the `one_plot` argument of [hill_pq()].
#' @param ... Additional arguments passed on to [patchwork::wrap_plots()]
#'   (e.g. `ncol`, `nrow`, `guides`) when `one_plot = TRUE`.
#'
#' @return If `one_plot = FALSE` (default), a named list of plot objects (the
#'   `alpha` entry is itself a multi-panel \pkg{patchwork} figure). If
#'   `one_plot = TRUE`, a single \code{\link[patchwork]{patchwork}} object
#'   assembling all panels.
#' @export
#' @author Adrien Taudière
#' @details The alpha-diversity panel relies on \pkg{patchwork} (always) and,
#'   for a numeric `fact`, on \pkg{ggstatsplot} (via [ggscatt_pq()]). Other
#'   panels may additionally require \pkg{ggVennDiagram} (Venn),
#'   \pkg{ComplexUpset} (UpSet) or \pkg{umap} (UMAP) depending on the options
#'   used.
#' @seealso [hill_bar_pq()], [ggscatt_pq()], [plot_ordination_pq()],
#'   [ggvenn_pq()], [upset_pq()], [umap_pq()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   # Balanced 6-sample subset spanning two Height levels (fast example).
#'   # UMAP and Venn are disabled here (too few samples / extra dependency).
#'   sn <- sample_names(data_fungi_mini)
#'   hi <- sn[which(data_fungi_mini@sam_data$Height == "High")[1:3]]
#'   lo <- sn[which(data_fungi_mini@sam_data$Height == "Low")[1:3]]
#'   ps <- prune_samples(c(hi, lo), data_fungi_mini)
#'   ps <- clean_pq(ps)
#'   plot_overview_pq(ps, fact = "Height", add_venn = FALSE, add_umap = FALSE)
#' }
#' }
#' \dontrun{
#' # Full overview with all panels (factor with 2 levels -> Venn).
#' plot_overview_pq(data_fungi_mini, fact = "Height", one_plot = TRUE)
#'
#' # Numeric variable -> Hill scatter (ggscatt_pq) + gradient-colored
#' # ordination/UMAP, Venn/UpSet automatically skipped.
#' plot_overview_pq(data_fungi_mini, fact = "Time", add_umap = FALSE)
#' }
plot_overview_pq <- function(
  physeq,
  fact,
  q = c(0, 1, 2),
  add_alpha = TRUE,
  add_ordination = TRUE,
  add_venn = TRUE,
  add_umap = TRUE,
  venn_max = 4,
  ordination_method = "PCoA",
  dist_method = "bray",
  one_plot = FALSE,
  ...
) {
  verify_pq(physeq)

  if (is.null(fact) || !fact %in% sample_variables(physeq)) {
    cli::cli_abort(
      "{.arg fact} must be a variable of the {.field sam_data} slot of {.arg physeq}."
    )
  }
  if (!add_alpha && !add_ordination && !add_venn && !add_umap) {
    cli::cli_abort(
      "At least one of {.arg add_alpha}, {.arg add_ordination}, {.arg add_venn} or {.arg add_umap} must be TRUE."
    )
  }
  # patchwork is needed to assemble the panels (one_plot = TRUE) and by the
  # alpha-diversity wrappers (hill_bar_pq / ggscatt_pq) for multi-panel output.
  if (
    (one_plot || add_alpha) && !requireNamespace("patchwork", quietly = TRUE)
  ) {
    cli::cli_abort(
      "Package {.pkg patchwork} is required by {.fn plot_overview_pq}. Please install it."
    )
  }

  # Drop samples with a missing modality so every panel sees the same samples.
  if (anyNA(physeq@sam_data[[fact]])) {
    physeq <- subset_samples_pq(physeq, !is.na(physeq@sam_data[[fact]]))
  }

  is_numeric_fact <- is.numeric(physeq@sam_data[[fact]])
  plots <- list()

  if (add_alpha) {
    if (is_numeric_fact) {
      plots[["alpha"]] <- ggscatt_pq(
        physeq,
        num_modality = fact,
        q = q,
        one_plot = TRUE,
        verbose = FALSE
      )
    } else {
      plots[["alpha"]] <- rlang::inject(
        hill_bar_pq(physeq, x = !!rlang::sym(fact), q = q)
      )
    }
  }

  if (add_venn) {
    if (is_numeric_fact) {
      cli::cli_inform(
        "{.arg fact} is numeric: the Venn/UpSet panel is skipped (it requires discrete groups)."
      )
    } else {
      n_lev <- nlevels(as.factor(physeq@sam_data[[fact]]))
      if (n_lev <= venn_max) {
        plots[["venn"]] <- ggvenn_pq(physeq, fact = fact, verbose = FALSE)
      } else {
        plots[["upset"]] <- upset_pq(physeq, fact = fact, verbose = FALSE)
      }
    }
  }

  if (add_ordination) {
    plots[["ordination"]] <- plot_ordination_pq(
      physeq,
      method = dist_method,
      ordination_method = ordination_method,
      color = fact
    )
  }

  if (add_umap) {
    # umap::umap defaults to n_neighbors = 15, which requires > 15 samples.
    if (!requireNamespace("umap", quietly = TRUE)) {
      cli::cli_inform(
        "The UMAP panel is skipped: package {.pkg umap} is not installed."
      )
    } else if (nsamples(physeq) <= 15) {
      cli::cli_inform(
        "The UMAP panel is skipped: {nsamples(physeq)} samples are too few for the default {.code n_neighbors = 15} of {.fn umap::umap}."
      )
    } else {
      df_umap <- umap_pq(physeq)
      plots[["umap"]] <- ggplot(
        df_umap,
        aes(
          x = .data[["x_umap"]],
          y = .data[["y_umap"]],
          color = .data[[fact]]
        )
      ) +
        geom_point(size = 2) +
        labs(title = "UMAP")
    }
  }

  if (one_plot) {
    return(patchwork::wrap_plots(plots, ...))
  }
  plots
}
################################################################################
