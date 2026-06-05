################################################################################
#' Fastp: adapter detection, quality control and preprocessing of FASTQ files
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' This function is a wrapper for the
#' [fastp](https://github.com/OpenGene/fastp) software to perform quality
#' control and preprocessing of FASTQ files (paired-end or single-end).
#'
#' To install fastp on Linux:
#'
#' ```bash
#' wget http://opengene.org/fastp/fastp
#' chmod a+x ./fastp
#' sudo mv fastp /usr/local/bin/
#' ```
#'
#' For other OS or specific installation instructions, please refer to the
#' official fastp GitHub repository: <https://github.com/OpenGene/fastp>.
#'
#' @param path_to_fastq (required) Path to the directory containing FASTQ files.
#'   See [list_fastq_files()] for help.
#' @inheritParams list_fastq_files
#' @param folder_output Name of the output folder (default: "fastp_output").
#' @param nproc Number of threads to use (default: 1).
#' @param paired_end Logical, whether files are paired-end (default: TRUE).
#' @param cmd_is_run Logical, whether to run the command (default: TRUE). If
#'   FALSE, the only effect of the function is to return a list of commands to
#'   run manually in a terminal.
#' @param return_file_path Logical, whether to return output file paths
#'   (default: FALSE). Useful in `targets` workflows.
#' @param qualified_quality_phred Minimum base quality (default: 20).
#' @param average_qual Minimum average read quality (default: 20).
#' @param length_required Minimum read length after trimming (default: 50).
#' @param detect_adapter Logical, enable automatic adapter detection for
#'   paired-end reads (default: TRUE).
#' @param trim_poly_g Logical, trim poly-G tails (default: TRUE).
#' @param correction Logical, correct mismatches in overlapped paired-end reads
#'   (default: TRUE).
#' @param cut_front Logical, trim low quality bases from 5' end (default: TRUE).
#' @param cut_tail Logical, trim low quality bases from 3' end (default: TRUE).
#' @param cut_window_size Window size for sliding quality filter (default: 4).
#' @param cut_mean_quality Minimum mean quality for trimming window
#'   (default: 20).
#' @param extra_fastp_args Additional fastp arguments as a character string
#'   (default: "").
#'
#' @return If `cmd_is_run = FALSE`, a named list of fastp commands. If
#'   `return_file_path = TRUE`, a list of output file paths. Otherwise `NULL`
#'   invisibly, after running fastp.
#' @export
#' @author Adrien Taudière
#' @seealso [is_fastp_installed()], [list_fastq_files()]
#' @details
#' This function is a wrapper for the fastp software, please cite fastp
#' (\doi{10.1093/bioinformatics/bty560}) if you use this function.
#'
#' @examples
#' \dontrun{
#' # Paired-end example
#' fastp(
#'   path_to_fastq = system.file("extdata", package = "MiscMetabar"),
#'   folder_output = file.path(tempdir(), "qc_data"),
#'   nproc = 8,
#'   paired_end = TRUE
#' )
#'
#' # Single-end example with extra options
#' fastp(
#'   path_to_fastq = "raw_data/",
#'   folder_output = file.path(tempdir(), "qc_data"),
#'   nproc = 8,
#'   paired_end = FALSE,
#'   extra_fastp_args = "--merge --n_base_limit 5 --dedup"
#' )
#'
#' # Only build the commands, do not run them
#' fastp(
#'   path_to_fastq = system.file("extdata", package = "MiscMetabar"),
#'   cmd_is_run = FALSE
#' )
#' }
fastp <- function(
  path_to_fastq,
  folder_output = "fastp_output",
  nproc = 1,
  pattern = "fastq.gz",
  pattern_R1 = "_R1",
  pattern_R2 = "_R2",
  nb_files = Inf,
  paired_end = TRUE,
  cmd_is_run = TRUE,
  return_file_path = FALSE,
  qualified_quality_phred = 20,
  average_qual = 20,
  length_required = 50,
  detect_adapter = TRUE,
  trim_poly_g = TRUE,
  correction = TRUE,
  cut_front = TRUE,
  cut_tail = TRUE,
  cut_window_size = 4,
  cut_mean_quality = 20,
  extra_fastp_args = ""
) {
  if (!dir.exists(folder_output)) {
    dir.create(folder_output, recursive = TRUE)
  }

  lff <- list_fastq_files(
    path_to_fastq,
    paired_end = paired_end,
    pattern = pattern,
    pattern_R1 = pattern_R1,
    pattern_R2 = pattern_R2,
    nb_files = nb_files
  )

  cmd <- vector("list", length(lff$fnfs)) # pre-allocated for performance
  names(cmd) <- lff$fnfs

  for (i in seq_along(lff$fnfs)) {
    f_R1 <- lff$fnfs[i]
    base_name <- sub("\\.fastq(\\.gz)?$", "", basename(f_R1))

    fastp_cmd <- paste0(
      "fastp",
      " -i ",
      f_R1,
      " -o ",
      folder_output,
      "/",
      basename(f_R1),
      " --thread ",
      nproc,
      " --qualified_quality_phred ",
      qualified_quality_phred,
      " --average_qual ",
      average_qual,
      " --length_required ",
      length_required,
      " --cut_window_size ",
      cut_window_size,
      " --cut_mean_quality ",
      cut_mean_quality
    )

    if (paired_end) {
      f_R2 <- lff$fnrs[i]
      fastp_cmd <- paste0(
        fastp_cmd,
        " -I ",
        f_R2,
        " -O ",
        folder_output,
        "/",
        basename(f_R2)
      )
      if (detect_adapter) {
        fastp_cmd <- paste0(fastp_cmd, " --detect_adapter_for_pe")
      }
      if (correction) {
        fastp_cmd <- paste0(fastp_cmd, " --correction")
      }
    }

    if (trim_poly_g) {
      fastp_cmd <- paste0(fastp_cmd, " --trim_poly_g")
    }
    if (cut_front) {
      fastp_cmd <- paste0(fastp_cmd, " --cut_front")
    }
    if (cut_tail) {
      fastp_cmd <- paste0(fastp_cmd, " --cut_tail")
    }

    fastp_cmd <- paste0(
      fastp_cmd,
      " --html ",
      folder_output,
      "/",
      base_name,
      ".fastp.html",
      " --json ",
      folder_output,
      "/",
      base_name,
      ".fastp.json"
    )

    if (nzchar(extra_fastp_args)) {
      fastp_cmd <- paste0(fastp_cmd, " ", extra_fastp_args)
    }

    cmd[[f_R1]] <- fastp_cmd
  }

  if (cmd_is_run) {
    script_path <- file.path(tempdir(), "script_fastp.sh")
    writeLines(unlist(cmd), script_path)
    system2("bash", script_path)
    message(
      "Output files are available in the folder ",
      normalizePath(folder_output)
    )
    unlink(script_path)
  } else {
    return(cmd)
  }

  if (return_file_path) {
    output_R1 <- file.path(folder_output, basename(lff$fnfs))
    if (paired_end) {
      output_R2 <- file.path(folder_output, basename(lff$fnrs))
      return(list("fnfs" = output_R1, "fnrs" = output_R2))
    } else {
      return(list("fnfs" = output_R1))
    }
  }

  invisible(NULL)
}
################################################################################

################################################################################
#' Test if fastp is installed.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage.
#'
#' @param path (default: fastp) Path to fastp.
#' @export
#' @return A logical that say if fastp is installed.
#' @seealso [fastp()]
#' @examples
#' MiscMetabar::is_fastp_installed()
#' @author Adrien Taudière
is_fastp_installed <- function(path = "fastp") {
  return(
    !inherits(
      try(
        system2(path, "--version", stdout = TRUE, stderr = TRUE),
        silent = TRUE
      ),
      "try-error"
    )
  )
}
################################################################################
