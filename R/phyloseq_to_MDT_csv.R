################################################################################
#' Export a phyloseq object to GBIF MDT template CSV/TSV files
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Write the OTU table, the taxonomy table and the sample data of a phyloseq
#' object to separate delimited text files, one per table, following the
#' GBIF [Metabarcoding Data Toolkit (MDT)](https://mdt.gbif.org/) template
#' layout. This is the plain-text (CSV/TSV) counterpart of
#' [phyloseq_to_MDT_excel()], which writes a single multi-sheet `.xlsx` file.
#'
#' The MDT template expects the OTU table with **OTU IDs in rows and sample IDs
#' in columns** (sequence read counts in the cells), and the sample and
#' taxonomy tables keyed by a leading `id` column. The reference sequences,
#' when available, are appended as a `DNA_sequence` column of the taxonomy file
#' (Darwin Core DNA-derived-data extension term).
#'
#' When `check_dwc = TRUE` (the default), a lightweight Darwin Core compliance
#' check warns about recommended sample-level terms that are missing from
#' `sample_data` (`decimalLatitude`, `decimalLongitude`, `eventDate`). This is
#' a non-blocking helper, not a full validation of the GBIF MDT template.
#'
#' @param physeq (required) a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param path (character, default `"."`) Directory where the files are
#'   written. Created (recursively) if it does not exist.
#' @param prefix (character, default `""`) Optional prefix prepended to each
#'   output file name (e.g. `"data_fungi_"`).
#' @param sep (character, default `","`) Field separator. Use `"\t"` to write
#'   tab-separated (TSV) files, the format favored by the GBIF MDT validator.
#'   The file extension follows `sep` (`.csv` for `","`, `.tsv` otherwise).
#' @param check_dwc (logical, default TRUE) If TRUE, warn about recommended
#'   Darwin Core sample-level terms missing from `sample_data`.
#'
#' @return Invisibly returns a named character vector of the written file paths
#'   (`OTU_table`, `Taxonomy`, `Samples`).
#' @export
#' @author Adrien Taudière
#' @details See the GBIF
#' [Metabarcoding Data Toolkit](https://mdt.gbif.org/) for the expected input
#' format. The MDT accepts both TSV and XLSX uploads.
#' @seealso [phyloseq_to_MDT_excel()]
#'
#' @examples
#' \donttest{
#' out_dir <- file.path(tempdir(), "mdt_csv")
#' files <- phyloseq_to_MDT_csv(clean_pq(data_fungi_mini), path = out_dir)
#' file.exists(files)
#' unlink(out_dir, recursive = TRUE)
#' }
#' \dontrun{
#' # Tab-separated output (MDT-favored TSV)
#' phyloseq_to_MDT_csv(data_fungi_mini, path = "mdt", sep = "\t")
#' }
phyloseq_to_MDT_csv <- function(
  physeq,
  path = ".",
  prefix = "",
  sep = ",",
  check_dwc = TRUE
) {
  verify_pq(physeq)

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  ext <- if (identical(sep, ",")) {
    ".csv"
  } else {
    ".tsv"
  }

  write_one <- function(df, name) {
    out_file <- file.path(path, paste0(prefix, name, ext))
    utils::write.table(
      df,
      file = out_file,
      sep = sep,
      row.names = FALSE,
      quote = TRUE,
      na = ""
    )
    out_file
  }

  # OTU table: OTU IDs in rows, sample IDs in columns (MDT orientation)
  if (taxa_are_rows(physeq)) {
    otu_table_df <- as.data.frame(as(physeq@otu_table, "matrix"))
  } else {
    otu_table_df <- as.data.frame(t(as(physeq@otu_table, "matrix")))
  }
  otu_table_df <- cbind(id = rownames(otu_table_df), otu_table_df)
  files <- c(OTU_table = write_one(otu_table_df, "OTU_table"))

  # Taxonomy table (+ refseq as DNA_sequence), keyed by OTU id
  if (!is.null(physeq@tax_table)) {
    tax_table_df <- as.data.frame(as(physeq@tax_table, "matrix"))
    if (!is.null(physeq@refseq)) {
      tax_table_df[["DNA_sequence"]] <- as.character(physeq@refseq)
    }
    tax_table_df <- cbind(id = rownames(tax_table_df), tax_table_df)
    files["Taxonomy"] <- write_one(tax_table_df, "Taxonomy")
  }

  # Sample data, keyed by sample id
  if (!is.null(physeq@sam_data)) {
    sample_data_df <- as.data.frame(as(physeq@sam_data, "data.frame"))
    sample_data_df <- cbind(id = rownames(sample_data_df), sample_data_df)
    files["Samples"] <- write_one(sample_data_df, "Samples")

    if (check_dwc) {
      recommended <- c("decimalLatitude", "decimalLongitude", "eventDate")
      missing_terms <- setdiff(recommended, colnames(sample_data_df))
      if (length(missing_terms) > 0) {
        cli::cli_warn(c(
          "!" = "Recommended Darwin Core sample terms missing from {.field sample_data}: {.val {missing_terms}}.",
          "i" = "Add them before GBIF MDT submission if available."
        ))
      }
    }
  }

  message("MDT template files written to ", normalizePath(path))
  invisible(files)
}
################################################################################
