################################################################################
#' Convert FASTQ file(s) to FASTA file(s)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Convert one or several FASTQ files (optionally gzip-compressed) to FASTA
#' format. For each input file, the `@` sign starting a FASTQ sequence
#' header is replaced by the `>` sign required by the FASTA format, and the
#' quality lines are dropped. Output files keep the input base name with the
#' extension changed to `.fasta`.
#'
#' Inspired by \url{https://github.com/ambuvjyn/baseq}.
#'
#' @param fastq_files (character, required) One or more paths to FASTQ
#'   files. The `.fastq`, `.fq`, `.fastq.gz` and `.fq.gz` extensions are
#'   recognized (case-insensitive).
#' @param output_folder (character, default: `tempdir()`) Folder where the
#'   FASTA files are written. Created if it does not exist.
#' @param force (logical, default: FALSE) If FALSE, existing output files
#'   are not overwritten and an error is thrown instead.
#'
#' @return A character vector of the paths to the output FASTA files
#'   (invisibly named by the input files).
#'
#' @author Adrien Taudière
#' @export
#'
#' @examples
#' \donttest{
#' fq <- system.file("extdata", "ex.fastq", package = "MiscMetabar")
#' out <- fastq_to_fasta(fq, output_folder = tempdir())
#' readLines(out, n = 2)
#' }
fastq_to_fasta <- function(
  fastq_files,
  output_folder = tempdir(),
  force = FALSE
) {
  if (!is.character(fastq_files) || length(fastq_files) == 0) {
    cli::cli_abort("{.arg fastq_files} must be a non-empty character vector.")
  }
  missing_files <- fastq_files[!file.exists(fastq_files)]
  if (length(missing_files) > 0) {
    cli::cli_abort(c(
      "Some {.arg fastq_files} do not exist:",
      "x" = "{.file {missing_files}}"
    ))
  }
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  out_paths <- vapply(
    fastq_files,
    function(fastq_file) {
      base_name <- basename(fastq_file)
      fasta_name <- sub(
        "\\.(fastq|fq)(\\.gz)?$",
        ".fasta",
        base_name,
        ignore.case = TRUE
      )
      if (identical(fasta_name, base_name)) {
        fasta_name <- paste0(base_name, ".fasta")
      }
      fasta_path <- file.path(output_folder, fasta_name)

      if (file.exists(fasta_path) && !force) {
        cli::cli_abort(c(
          "Output file already exists: {.file {fasta_path}}.",
          "i" = "Use {.code force = TRUE} to overwrite."
        ))
      }

      con <- if (grepl("\\.gz$", fastq_file, ignore.case = TRUE)) {
        gzfile(fastq_file, open = "rt")
      } else {
        file(fastq_file, open = "rt")
      }
      on.exit(close(con), add = TRUE)
      fastq_data <- readLines(con)

      if (length(fastq_data) %% 4 != 0) {
        cli::cli_warn(c(
          "{.file {fastq_file}} does not have a number of lines divisible by 4.",
          "i" = "It may not be a valid FASTQ file; output may be incomplete."
        ))
      }

      idx <- seq(1, length(fastq_data), by = 4)
      headers <- sub("^@", ">", fastq_data[idx])
      seqs <- gsub(" ", "", fastq_data[idx + 1])
      writeLines(c(rbind(headers, seqs)), fasta_path)

      fasta_path
    },
    FUN.VALUE = character(1)
  )

  cli::cli_inform(c(
    "v" = "Wrote {length(out_paths)} FASTA file{?s} to {.file {output_folder}}."
  ))
  invisible(out_paths)
}
