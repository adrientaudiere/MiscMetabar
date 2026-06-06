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
#' The function is column-name agnostic: point `fc` and `padj` at the relevant
#' columns of your result table (see *Details* for common cases).
#'
#' @param df (required) A `data.frame` (or any object coercible with
#'   [as.data.frame()], such as a `DESeqResults` object) with one row per taxon
#'   and at least a fold-change column and an adjusted p-value column.
#' @param fc (character, default `"log2FoldChange"`) Name of the (log2) fold
#'   change column.
#' @param padj (character, default `"padj"`) Name of the adjusted p-value
#'   column.
#' @param alpha (numeric, default 0.05) Adjusted p-value threshold for
#'   significance (horizontal dashed line).
#' @param lfc_threshold (numeric, default 1) Absolute (log2) fold-change
#'   threshold for biological relevance (vertical dashed lines). Set to 0 to
#'   classify on significance only.
#' @param label_col (character, default NULL) Optional column used to label the
#'   significant points (e.g. a taxonomic rank). If NULL, no labels are drawn.
#' @param label_n (integer, default 10) Maximum number of significant points to
#'   label, ranked by \eqn{-\log_{10}(\text{padj})}. Ignored when `label_col`
#'   is NULL.
#' @param point_size (numeric, default 2) Size of the points.
#' @param point_alpha (numeric, default 0.7) Opacity of the points.
#' @param palette (named character vector) Colours for the three statuses. Must
#'   be named `Down`, `NotDA` and `Up`.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
#' @author Adrien Taudière
#' @details
#' Common column mappings:
#' - \pkg{DESeq2} ([DESeq2::results()]): `fc = "log2FoldChange"`,
#'   `padj = "padj"` (the defaults).
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
#' # From a real DESeq2 analysis
#' data("GlobalPatterns", package = "phyloseq")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#' GP <- subset_samples(GP, SampleType %in% c("Soil", "Skin"))
#' if (requireNamespace("DESeq2")) {
#'   dds <- DESeq2::DESeq(phyloseq_to_deseq2(GP, ~SampleType))
#'   plot_volcano_pq(DESeq2::results(dds))
#' }
#' }
plot_volcano_pq <- function(
  df,
  fc = "log2FoldChange",
  padj = "padj",
  alpha = 0.05,
  lfc_threshold = 1,
  label_col = NULL,
  label_n = 10,
  point_size = 2,
  point_alpha = 0.7,
  palette = c(Down = "#1d2949", NotDA = "grey70", Up = "#4B3E1E")
) {
  df <- as.data.frame(df)

  if (!fc %in% colnames(df)) {
    cli::cli_abort(
      "Fold-change column {.val {fc}} not found in {.arg df}. Available columns: {.val {colnames(df)}}."
    )
  }
  if (!padj %in% colnames(df)) {
    cli::cli_abort(
      "Adjusted p-value column {.val {padj}} not found in {.arg df}. Available columns: {.val {colnames(df)}}."
    )
  }

  fc_vals <- as.numeric(df[[fc]])
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
