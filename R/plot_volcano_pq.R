################################################################################
#' Volcano plot for differential abundance results
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Build a volcano plot (log2 fold change on the x-axis against
#' \eqn{-\log_{10}(\text{adjusted p-value})} on the y-axis) from the result
#' table of a differential abundance analysis, e.g. \pkg{DESeq2}
#' ([DESeq2::results()]), \pkg{ANCOMBC} ([ancombc_pq()]) or \pkg{ALDEx2}
#' ([aldex_pq()]). Each taxon is colored according to its status: up, down or
#' not differentially abundant (`NotDA`).
#'
#' The function auto-detects the input type and sets `fc` and `padj`
#' accordingly (see *Details*). Pass explicit values to override.
#'
#' For \pkg{lefser} results ([lefser_pq()]), which have no per-feature
#' p-value, a score plot is produced instead: LDA score on both axes
#' (\eqn{|\text{score}|} on y), with taxa classified by sign and
#' `lfc_threshold`.
#'
#' @param df (required) A `data.frame` (or any object coercible with
#'   [as.data.frame()], such as a `DESeqResults` object), the `$res` list
#'   returned by [ancombc_pq()], or a `lefser_df` object from [lefser_pq()].
#' @param fc (character or NULL, default `NULL`) Name of the (log2) fold
#'   change column. When `NULL` the column is inferred from the input type
#'   (see *Details*).
#' @param padj (character or NULL, default `NULL`) Name of the adjusted
#'   p-value column. When `NULL` the column is inferred from the input type.
#'   Set to `NULL` explicitly (and supply a `lefser_df`) to activate score
#'   mode.
#' @param alpha (numeric, default 0.05) Adjusted p-value threshold for
#'   significance (horizontal dashed line). Ignored in LEfSe score mode.
#' @param lfc_threshold (numeric, default 1) Absolute (log2) fold-change
#'   threshold for biological relevance (vertical dashed lines). Set to 0 to
#'   classify on significance only.
#' @param label_col (character, default NULL) Optional column used to label the
#'   significant points (e.g. a taxonomic rank). If NULL, no labels are drawn.
#' @param label_n (integer, default 10) Maximum number of significant points to
#'   label, ranked by \eqn{-\log_{10}(\text{padj})} (or \eqn{|\text{score}|}
#'   in LEfSe mode). Ignored when `label_col` is NULL.
#' @param point_size (numeric, default 2) Size of the points.
#' @param point_alpha (numeric, default 0.7) Opacity of the points.
#' @param palette (named character vector) Colours for the three statuses. Must
#'   be named `Down`, `NotDA` and `Up`.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
#' @author Adrien Taudière
#' @details
#' **Auto-detection rules** (applied when `fc` or `padj` is `NULL`):
#'
#' - If `df` is a list with a `$res` slot containing `lfc_*` and `q_*`
#'   columns (output of [ancombc_pq()]): `$res` is extracted automatically,
#'   and the first `lfc_*` / `q_*` column is used. Pass `fc` / `padj`
#'   explicitly to pick a specific comparison when multiple groups are present.
#' - If `df` inherits from `lefser_df` (output of [lefser_pq()]): score mode
#'   is activated (`fc = "scores"`, no p-value axis).
#' - If columns `log2FoldChange` and `padj` are present (e.g. after
#'   `as.data.frame(DESeq2::results(dds))`): DESeq2 defaults are used.
#' - If columns `effect` and `wi.eBH` are present (output of [aldex_pq()]):
#'   ALDEx2 defaults are used.
#'
#' **Manual column mappings** (for other tools or custom results):
#' - \pkg{DESeq2} ([DESeq2::results()]): `fc = "log2FoldChange"`,
#'   `padj = "padj"`.
#' - \pkg{ANCOMBC} ([ancombc_pq()] `$res`): `fc = "lfc_<group>"`,
#'   `padj = "q_<group>"`.
#' - \pkg{ALDEx2} ([aldex_pq()]): `fc = "effect"`, `padj = "wi.eBH"` (or
#'   `"we.eBH"`).
#'
#' Taxa with `padj == 0` are drawn at the top of the plot (their
#' \eqn{-\log_{10}} value is capped just above the largest finite value); taxa
#' with a missing `padj` are classified as `NotDA`.
#'
#' @examples
#' \donttest{
#' # Synthetic differential abundance table (DESeq2-like columns)
#' set.seed(42)
#' res <- data.frame(
#'   log2FoldChange = rnorm(200, sd = 2),
#'   padj = runif(200)^3,
#'   Genus = sample(paste0("Genus_", 1:20), 200, replace = TRUE)
#' )
#' plot_volcano_pq(res)
#' plot_volcano_pq(res, lfc_threshold = 2, label_col = "Genus", label_n = 5)
#' }
#' \dontrun{
#' # From a real DESeq2 analysis (auto-detected)
#' data("GlobalPatterns", package = "phyloseq")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#' GP <- subset_samples(GP, SampleType %in% c("Soil", "Skin"))
#' if (requireNamespace("DESeq2")) {
#'   dds <- DESeq2::DESeq(phyloseq_to_deseq2(GP, ~SampleType))
#'   plot_volcano_pq(DESeq2::results(dds))
#' }
#'
#' # From ANCOMBC (auto-detected: $res extracted, lfc_/q_ columns picked)
#' if (requireNamespace("ANCOMBC")) {
#'   res_ancombc <- ancombc_pq(data_fungi_mini, fact=  "Height",  levels_fact = c("Low", "High"))
#'   plot_volcano_pq(res_ancombc)
#' }
#'
#' # From ALDEx2 (auto-detected)
#' if (requireNamespace("ALDEx2")) {
#'   res_aldex <- aldex_pq(data_fungi_mini, bifactor = "Height",
#'                         modalities = c("Low", "High"))
#'   plot_volcano_pq(res_aldex)
#' }
#'
#' # From LEfSe (score plot, no p-value)
#' if (requireNamespace("lefser") && requireNamespace("mia")) {
#'   res_lefse <- lefser_pq(data_fungi_mini, bifactor = "Height",
#'                          modalities = c("Low", "High"))
#'   plot_volcano_pq(res_lefse)
#' }
#' }
plot_volcano_pq <- function(
  df,
  fc = NULL,
  padj = NULL,
  alpha = 0.05,
  lfc_threshold = 1,
  label_col = NULL,
  label_n = 10,
  point_size = 2,
  point_alpha = 0.7,
  palette = c(Down = "#1d2949", NotDA = "grey70", Up = "#4B3E1E")
) {
  # Auto-detect ANCOMBC: raw list with $res slot containing lfc_/q_ columns
  if (
    is.list(df) && !is.data.frame(df) &&
      "res" %in% names(df) && is.data.frame(df$res) &&
      any(grepl("^lfc_", colnames(df$res)))
  ) {
    df <- df$res
    if (is.null(fc)) {
      lfc_cols <- grep("^lfc_", colnames(df), value = TRUE)
      lfc_cols <- lfc_cols[!grepl("Intercept", lfc_cols, fixed = TRUE)]
      if (length(lfc_cols) > 1) {
        cli::cli_inform(
          "Multiple lfc columns detected ({.val {lfc_cols}}); using {.val {lfc_cols[1]}}. Set {.arg fc} to override."
        )
      }
      if (length(lfc_cols) >= 1) {
        fc <- lfc_cols[1]
      }
    }
    if (is.null(padj)) {
      q_cols <- grep("^q_", colnames(df), value = TRUE)
      q_cols <- q_cols[!grepl("Intercept", q_cols, fixed = TRUE)]
      if (length(q_cols) > 1) {
        cli::cli_inform(
          "Multiple q columns detected ({.val {q_cols}}); using {.val {q_cols[1]}}. Set {.arg padj} to override."
        )
      }
      if (length(q_cols) >= 1) {
        padj <- q_cols[1]
      }
    }
  }

  # Auto-detect LEfSe: lefser_df class has scores but no p-value
  lefser_mode <- inherits(df, "lefser_df")
  if (lefser_mode) {
    if (is.null(fc)) {
      fc <- "scores"
    }
    padj <- NULL
  }

  df <- as.data.frame(df)

  # Auto-detect DESeq2 columns (after coercion, DESeqResults → data.frame)
  if (is.null(fc) && "log2FoldChange" %in% colnames(df)) {
    fc <- "log2FoldChange"
  }
  if (is.null(padj) && !lefser_mode && "padj" %in% colnames(df)) {
    padj <- "padj"
  }

  # Auto-detect ALDEx2 columns
  if (is.null(fc) && "effect" %in% colnames(df)) {
    fc <- "effect"
  }
  if (is.null(padj) && !lefser_mode && "wi.eBH" %in% colnames(df)) {
    padj <- "wi.eBH"
  }

  # Ultimate fallback (preserves old default behaviour)
  if (is.null(fc)) {
    fc <- "log2FoldChange"
  }
  if (is.null(padj) && !lefser_mode) {
    padj <- "padj"
  }

  if (!fc %in% colnames(df)) {
    cli::cli_abort(
      "Fold-change column {.val {fc}} not found in {.arg df}. Available columns: {.val {colnames(df)}}."
    )
  }
  if (!is.null(padj) && !padj %in% colnames(df)) {
    cli::cli_abort(
      "Adjusted p-value column {.val {padj}} not found in {.arg df}. Available columns: {.val {colnames(df)}}."
    )
  }

  fc_vals <- as.numeric(df[[fc]])

  if (lefser_mode) {
    # Score plot for LEfSe: classify by sign × threshold, y = |score|
    status <- rep("NotDA", nrow(df))
    status[fc_vals >= lfc_threshold] <- "Up"
    status[fc_vals <= -lfc_threshold] <- "Down"
    status <- factor(status, levels = c("Down", "NotDA", "Up"))

    plot_df <- data.frame(
      .fc = fc_vals,
      .neglog10 = abs(fc_vals),
      .status = status,
      stringsAsFactors = FALSE
    )
    if (!is.null(label_col)) {
      if (!label_col %in% colnames(df)) {
        cli::cli_abort("Label column {.val {label_col}} not found in {.arg df}.")
      }
      plot_df$.label <- as.character(df[[label_col]])
    }

    p <- ggplot(
      plot_df,
      aes(x = .data$.fc, y = .data$.neglog10, color = .data$.status)
    ) +
      geom_vline(
        xintercept = c(-lfc_threshold, lfc_threshold),
        linetype = "dashed",
        color = "grey60"
      ) +
      geom_point(size = point_size, alpha = point_alpha) +
      scale_color_manual(values = palette, drop = FALSE, name = NULL) +
      labs(x = "LDA score", y = "|LDA score|") +
      theme_bw()

    if (!is.null(label_col)) {
      lab_df <- plot_df[plot_df$.status != "NotDA", , drop = FALSE]
      if (nrow(lab_df) > 0) {
        lab_df <- lab_df[order(-lab_df$.neglog10), , drop = FALSE]
        lab_df <- utils::head(lab_df, label_n)
        p <- p +
          geom_text(
            data = lab_df,
            aes(label = .data$.label),
            size = 3,
            vjust = -0.6,
            show.legend = FALSE
          )
      }
    }
    return(p)
  }

  # Standard volcano mode (DESeq2, ALDEx2, ANCOMBC, ...)
  padj_vals <- as.numeric(df[[padj]])

  # Significance status
  is_signif <- !is.na(padj_vals) & padj_vals < alpha
  status <- rep("NotDA", nrow(df))
  status[is_signif & fc_vals >= lfc_threshold] <- "Up"
  status[is_signif & fc_vals <= -lfc_threshold] <- "Down"
  status <- factor(status, levels = c("Down", "NotDA", "Up"))

  # -log10(padj), capping padj == 0 just above the largest finite value
  neglog10 <- -log10(padj_vals)
  finite_max <- suppressWarnings(max(neglog10[is.finite(neglog10)]))
  if (any(is.infinite(neglog10)) && is.finite(finite_max)) {
    cli::cli_inform(
      "Some {.arg {padj}} values equal 0; capping their -log10 at {round(finite_max * 1.05, 2)}."
    )
    neglog10[is.infinite(neglog10)] <- finite_max * 1.05
  }
  neglog10[is.na(neglog10)] <- 0

  plot_df <- data.frame(
    .fc = fc_vals,
    .neglog10 = neglog10,
    .status = status,
    stringsAsFactors = FALSE
  )
  if (!is.null(label_col)) {
    if (!label_col %in% colnames(df)) {
      cli::cli_abort(
        "Label column {.val {label_col}} not found in {.arg df}."
      )
    }
    plot_df$.label <- as.character(df[[label_col]])
  }

  p <- ggplot(
    plot_df,
    aes(x = .data$.fc, y = .data$.neglog10, color = .data$.status)
  ) +
    geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = palette, drop = FALSE, name = NULL) +
    labs(
      x = expression(log[2] ~ "fold change"),
      y = expression(-log[10] ~ "(adjusted p-value)")
    ) +
    theme_bw()

  # Optional labels on the most significant points
  if (!is.null(label_col)) {
    lab_df <- plot_df[plot_df$.status != "NotDA", , drop = FALSE]
    if (nrow(lab_df) > 0) {
      lab_df <- lab_df[order(-lab_df$.neglog10), , drop = FALSE]
      lab_df <- utils::head(lab_df, label_n)
      p <- p +
        geom_text(
          data = lab_df,
          aes(label = .data$.label),
          size = 3,
          vjust = -0.6,
          show.legend = FALSE
        )
    }
  }

  return(p)
}
################################################################################
