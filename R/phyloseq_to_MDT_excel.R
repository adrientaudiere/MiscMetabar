################################################################################
#' Export a phyloseq object to a multi-sheet Excel file for GBIF MDT submission
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Write the OTU table, the sample data, the taxonomy table and (if present) the
#' reference sequences of a phyloseq object to a single multi-sheet `.xlsx`
#' file, formatted for submission to the GBIF
#' [Metabarcoding Data Toolkit (MDT)](https://mdt.gbif.org/).
#'
#' Each sheet gets a leading `id` column holding the row identifiers (taxa or
#' samples names). The reference sequences, when available, are appended as a
#' `DNA_sequence` column of the taxonomy sheet (Darwin Core DNA-derived-data
#' extension term).
#'
#' When `check_dwc = TRUE` (the default), a lightweight Darwin Core compliance
#' check warns about recommended sample-level terms that are missing from
#' `sample_data` (`decimalLatitude`, `decimalLongitude`, `eventDate`). This is
#' a non-blocking helper, not a full validation of the GBIF MDT template.
#'
#' @param physeq (required) a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param filename (character, default `"Phyloseq_Tables.xlsx"`) Path of the
#'   output `.xlsx` file.
#' @param check_dwc (logical, default TRUE) If TRUE, warn about recommended
#'   Darwin Core sample-level terms missing from `sample_data`.
#'
#' @return Invisibly returns the path to the written file (`filename`).
#' @export
#' @author Adrien Taudière
#' @details
#' This function requires the \pkg{writexl} package. See the GBIF
#' [Metabarcoding Data Toolkit](https://mdt.gbif.org/) for the expected
#' input format.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("writexl")) {
#'   out <- file.path(tempdir(), "data_fungi_mini_MDT.xlsx")
#'   phyloseq_to_MDT_excel(clean_pq(data_fungi_mini), filename = out)
#'   file.exists(out)
#'   unlink(out)
#' }
#' }
phyloseq_to_MDT_excel <- function(
  physeq,
  filename = "Phyloseq_Tables.xlsx",
  check_dwc = TRUE
) {
  if (!requireNamespace("writexl", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg writexl} is required for {.fn phyloseq_to_MDT_excel}. Please install it."
    )
  }
  verify_pq(physeq)

  # OTU table: samples in rows, taxa in columns
  if (taxa_are_rows(physeq)) {
    otu_table_df <- as.data.frame(t(as(physeq@otu_table, "matrix")))
  } else {
    otu_table_df <- as.data.frame(as(physeq@otu_table, "matrix"))
  }
  otu_sheet <- cbind(id = rownames(otu_table_df), otu_table_df)

  # Sample data
  sample_sheet <- NULL
  if (!is.null(physeq@sam_data)) {
    sample_data_df <- as.data.frame(as(physeq@sam_data, "data.frame"))
    sample_sheet <- cbind(id = rownames(sample_data_df), sample_data_df)

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

  # Taxonomy table (+ refseq as DNA_sequence)
  taxonomy_sheet <- NULL
  if (!is.null(physeq@tax_table)) {
    tax_table_df <- as.data.frame(as(physeq@tax_table, "matrix"))
    if (!is.null(physeq@refseq)) {
      tax_table_df[["DNA_sequence"]] <- as.character(physeq@refseq)
    }
    taxonomy_sheet <- cbind(id = rownames(tax_table_df), tax_table_df)
  }

  sheets <- list(OTU_table = otu_sheet)
  if (!is.null(sample_sheet)) {
    sheets[["Samples"]] <- sample_sheet
  }
  if (!is.null(taxonomy_sheet)) {
    sheets[["Taxonomy"]] <- taxonomy_sheet
  }

  writexl::write_xlsx(sheets, path = filename)
  message("Excel file written to ", normalizePath(filename))
  invisible(filename)
}
################################################################################
