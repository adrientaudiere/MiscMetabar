################################################################################
#' Find the ITSxpress executable
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Resolution order: the `MiscMetabar.itsxpresspath` option (if set), then the
#' system `PATH`. ITSxpress is usually installed in a conda environment (see
#' [install_itsxpress()]); in that case it is not on the `PATH` of R and the
#' environment is activated with the `args_before_itsxpress` argument of
#' [itsxpress_pq()] and [is_itsxpress_installed()] instead.
#'
#' @return A character string with the path to `itsxpress`, or `"itsxpress"` as
#'   a fallback (relying on `PATH` resolution, possibly after
#'   `args_before_itsxpress`).
#' @export
#' @examples
#' find_itsxpress()
#' @author Adrien Taudière
#' @seealso [itsxpress_pq()], [is_itsxpress_installed()], [install_itsxpress()]
find_itsxpress <- function() {
  opt <- getOption("MiscMetabar.itsxpresspath")
  if (!is.null(opt) && nzchar(opt)) {
    return(opt)
  }
  on_path <- unname(Sys.which("itsxpress"))
  if (nzchar(on_path)) {
    return(on_path)
  }
  "itsxpress"
}

################################################################################
#' Test if ITSxpress is installed
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful for `testthat` and examples compilation for `R CMD check` and
#' test coverage.
#'
#' @param path (default: [find_itsxpress()]) Path to `itsxpress`.
#' @param args_before_itsxpress (String, default `""`) A one-line bash command
#'   run before `itsxpress`, e.g. the value returned by [install_itsxpress()]:
#'   `"source ~/miniforge3/etc/profile.d/conda.sh && conda activate itsxpressenv && "`.
#' @export
#' @return A logical indicating whether `itsxpress` can be run.
#'
#' @examples
#' MiscMetabar::is_itsxpress_installed()
#' @author Adrien Taudière
#' @seealso [find_itsxpress()], [install_itsxpress()], [itsxpress_pq()]
is_itsxpress_installed <- function(
  path = find_itsxpress(),
  args_before_itsxpress = ""
) {
  status <- suppressWarnings(system2(
    "bash",
    c("-c", shQuote(paste0(args_before_itsxpress, shQuote(path), " --help"))),
    stdout = FALSE,
    stderr = FALSE
  ))
  identical(as.integer(status), 0L)
}

################################################################################
#' Install ITSxpress in a conda environment
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' ITSxpress (Rivers et al. 2018; Einarsson & Rivers 2024) is a Python tool
#' that requires VSEARCH, HMMER 3 and Biopython. The simplest way to install
#' all of them is the bioconda package, which this function installs in a
#' dedicated conda environment. It returns the string to pass as
#' `args_before_itsxpress` to [itsxpress_pq()] and [is_itsxpress_installed()],
#' which activates that environment before each call.
#'
#' @param env_name (default: `"itsxpressenv"`) Name of the conda environment.
#' @param conda (default: `NULL`) Path to `mamba` or `conda`. If `NULL`,
#'   `mamba` then `conda` are looked up on the `PATH`.
#' @param force (default: `FALSE`) If `TRUE`, the environment is created even
#'   if one of that name already exists (conda then replaces it).
#'
#' @return The `args_before_itsxpress` string activating the environment
#'   (invisibly).
#'
#' @export
#' @examples
#' \dontrun{
#' prelude <- install_itsxpress()
#' is_itsxpress_installed(args_before_itsxpress = prelude)
#' }
#' @author Adrien Taudière
#' @seealso [itsxpress_pq()], [is_itsxpress_installed()], [find_itsxpress()]
install_itsxpress <- function(
  env_name = "itsxpressenv",
  conda = NULL,
  force = FALSE
) {
  if (is.null(conda)) {
    candidates <- unname(Sys.which(c("mamba", "conda")))
    conda <- candidates[nzchar(candidates)][1]
  }
  if (is.na(conda) || !nzchar(conda)) {
    cli::cli_abort(c(
      "Neither {.code mamba} nor {.code conda} was found on the PATH.",
      "i" = "Install Miniforge (\\url{https://github.com/conda-forge/miniforge})
      or install {.code itsxpress} via bioconda/pip and set
      {.code options(MiscMetabar.itsxpresspath = '/path/to/itsxpress')}."
    ))
  }
  conda_base <- system2(conda, c("info", "--base"), stdout = TRUE)
  conda_base <- conda_base[length(conda_base)]
  prelude <- paste0(
    "source ",
    shQuote(file.path(conda_base, "etc", "profile.d", "conda.sh")),
    " && conda activate ",
    env_name,
    " && "
  )

  if (!force && is_itsxpress_installed(args_before_itsxpress = prelude)) {
    message(
      "ITSxpress is already installed in the conda environment '",
      env_name,
      "'. Use force = TRUE to reinstall."
    )
    return(invisible(prelude))
  }

  message("Creating the conda environment '", env_name, "' with ITSxpress...")
  status <- system2(
    conda,
    c(
      "create",
      "-y",
      "-n",
      env_name,
      "-c",
      "conda-forge",
      "-c",
      "bioconda",
      "itsxpress"
    )
  )
  if (status != 0 || !is_itsxpress_installed(args_before_itsxpress = prelude)) {
    cli::cli_abort(
      "The installation of ITSxpress in the conda environment {.val {env_name}} failed."
    )
  }
  message("ITSxpress installed. Use args_before_itsxpress = \"", prelude, "\"")
  invisible(prelude)
}

# Allowed taxonomic groups for ITSxpress --taxa
itsxpress_taxa_groups <- c(
  "Alveolata",
  "Bryophyta",
  "Bacillariophyta",
  "Amoebozoa",
  "Euglenozoa",
  "Fungi",
  "Chlorophyta",
  "Rhodophyta",
  "Phaeophyceae",
  "Marchantiophyta",
  "Metazoa",
  "Oomycota",
  "Haptophyceae",
  "Raphidophyceae",
  "Rhizaria",
  "Synurophyceae",
  "Tracheophyta",
  "Eustigmatophyceae",
  "Parabasalia",
  "All"
)

################################################################################
#' Trim FASTQ files to the ITS region with ITSxpress before DADA2 denoising
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Run [ITSxpress](https://github.com/USDA-ARS-GBRU/itsxpress) (Rivers et al.
#' 2018; Einarsson & Rivers 2024) on paired-end or single-end FASTQ files to
#' trim off the conserved ribosomal flanks (SSU, 5.8S, LSU) while preserving
#' per-base Phred quality scores. Unlike [itsx_pq()], which runs ITSx on
#' already-denoised FASTA reference sequences in a `phyloseq` object,
#' `itsxpress_pq()` operates directly on raw FASTQ reads **before** DADA2
#' denoising ([filter_trim()], [denoised_reads()]).
#'
#' By default (`paired_end = TRUE` and `single_file_output = FALSE`),
#' `itsxpress_pq()` writes unmerged forward and reverse trimmed FASTQ files
#' (`--outfile` and `--outfile2`), which is the layout required for DADA2's
#' paired-end error model.
#'
#' @param path_to_fastq (required) Path to the directory containing FASTQ
#'   files. See [list_fastq_files()] for help.
#' @inheritParams list_fastq_files
#' @param folder_output (default: `"itsxpress_output"`) Path to the output
#'   directory where trimmed FASTQ files and per-sample logs are written.
#' @param region (default: `"ITS1"`) The ITS region to keep: `"ITS1"`,
#'   `"ITS2"` or `"ALL"` (full ITS1 + 5.8S + ITS2).
#' @param taxa (default: `"Fungi"`) Taxonomic group sequenced (`--taxa`). One
#'   of `"Fungi"`, `"All"`, `"Alveolata"`, `"Bryophyta"`, `"Bacillariophyta"`,
#'   `"Amoebozoa"`, `"Euglenozoa"`, `"Chlorophyta"`, `"Rhodophyta"`,
#'   `"Phaeophyceae"`, `"Marchantiophyta"`, `"Metazoa"`, `"Oomycota"`,
#'   `"Haptophyceae"`, `"Raphidophyceae"`, `"Rhizaria"`, `"Synurophyceae"`,
#'   `"Tracheophyta"`, `"Eustigmatophyceae"`, or `"Parabasalia"`.
#' @param paired_end (logical, default `TRUE`) Whether the input FASTQ files
#'   are paired-end.
#' @param single_file_output (logical, default `FALSE`) Only used when
#'   `paired_end = TRUE`. If `FALSE` (default, required for downstream DADA2
#'   paired-end denoising), `itsxpress` outputs unmerged R1 and R2 files via
#'   `--outfile` and `--outfile2`. If `TRUE`, `itsxpress` outputs a single
#'   merged FASTQ file per pair via `--outfile`.
#' @param cluster_id (numeric, default `1.0`) Percent identity in `[0.99, 1.0]`
#'   for clustering reads with VSEARCH (`--cluster_id`). Set to `1.0` (default)
#'   for exact dereplication when calling ASVs.
#' @param reversed_primers (logical, default `FALSE`) If `TRUE`, pass
#'   `--reversed_primers` for primers in reverse orientation (Taylor et al.
#'   2016).
#' @param trim_ccs (logical, default `FALSE`) If `TRUE`, pass `--trim-ccs` for
#'   PacBio CCS long reads (only valid when `paired_end = FALSE`).
#' @param nproc (integer, default `1`) Number of processor threads (`--threads`).
#' @param cmd_is_run (logical, default `TRUE`) Whether to execute the
#'   `itsxpress` commands. If `FALSE`, the function only builds and returns a
#'   named list of shell commands.
#' @param return_file_path (logical, default `FALSE`) If `TRUE` (when
#'   `cmd_is_run = TRUE`), return a list of output FASTQ file paths (`fnfs`
#'   and optionally `fnrs`), ready to be chained in a `targets` pipeline.
#' @param itsxpresspath (default: [find_itsxpress()]) Path to the `itsxpress`
#'   executable.
#' @param args_before_itsxpress (String, default `""`) A one-line bash command
#'   run before `itsxpress`, e.g. the conda activation string returned by
#'   [install_itsxpress()].
#' @param extra_itsxpress_args (String, default `""`) Additional command-line
#'   flags passed verbatim to `itsxpress`.
#'
#' @return If `cmd_is_run = FALSE`, a named list of `itsxpress` shell commands.
#'   If `return_file_path = TRUE`, a named list of output file paths (`fnfs`,
#'   and `fnrs` when `paired_end = TRUE` and `single_file_output = FALSE`).
#'   Otherwise `NULL` invisibly, after running `itsxpress`.
#' @export
#' @author Adrien Taudière
#' @references
#' Rivers AR, Weber KC, Gardner TG, Liu S, Armstrong SD (2018). ITSxpress:
#' Software to rapidly trim internally transcribed spacer sequences with
#' quality scores for marker gene analysis. *F1000Research* 7: 1418.
#' \doi{10.12688/f1000research.15704.1}
#'
#' Einarsson SV, Rivers AR (2024). ITSxpress version 2: software to rapidly
#' trim internal transcribed spacer sequences with quality scores for amplicon
#' sequencing. *Microbiology Spectrum* 12(4): e00601-24.
#' \doi{10.1128/spectrum.00601-24}
#' @seealso [install_itsxpress()], [is_itsxpress_installed()],
#'   [find_itsxpress()], [itsx_pq()], [fastp()], [cutadapt_remove_primers()]
#' @examples
#' # Build paired-end ITSxpress commands without running them
#' itsxpress_pq(
#'   path_to_fastq = system.file("extdata", package = "MiscMetabar"),
#'   folder_output = file.path(tempdir(), "itsxpress_out"),
#'   region = "ITS2",
#'   cmd_is_run = FALSE
#' )
#' \dontrun{
#' prelude <- install_itsxpress()
#' itsxpress_pq(
#'   path_to_fastq = system.file("extdata", package = "MiscMetabar"),
#'   folder_output = file.path(tempdir(), "itsxpress_out"),
#'   region = "ITS1",
#'   taxa = "Fungi",
#'   nproc = 2,
#'   args_before_itsxpress = prelude
#' )
#' }
itsxpress_pq <- function(
  path_to_fastq,
  folder_output = "itsxpress_output",
  region = c("ITS1", "ITS2", "ALL"),
  taxa = "Fungi",
  paired_end = TRUE,
  single_file_output = FALSE,
  cluster_id = 1.0,
  reversed_primers = FALSE,
  trim_ccs = FALSE,
  nproc = 1,
  pattern = "fastq.gz",
  pattern_R1 = "_R1",
  pattern_R2 = "_R2",
  nb_files = Inf,
  cmd_is_run = TRUE,
  return_file_path = FALSE,
  itsxpresspath = find_itsxpress(),
  args_before_itsxpress = "",
  extra_itsxpress_args = ""
) {
  region <- match.arg(region)
  taxa <- match.arg(taxa, itsxpress_taxa_groups)
  if (
    !is.numeric(cluster_id) ||
      length(cluster_id) != 1 ||
      is.na(cluster_id) ||
      cluster_id < 0.99 ||
      cluster_id > 1.0
  ) {
    cli::cli_abort(
      "{.arg cluster_id} must be a single number between 0.99 and 1.0, not {.val {cluster_id}}."
    )
  }
  if (trim_ccs && paired_end) {
    cli::cli_abort(
      "{.arg trim_ccs} = TRUE can only be used with single-end reads ({.code paired_end = FALSE})."
    )
  }

  lff <- list_fastq_files(
    path_to_fastq,
    paired_end = paired_end,
    pattern = pattern,
    pattern_R1 = pattern_R1,
    pattern_R2 = pattern_R2,
    nb_files = nb_files
  )

  cmd <- vector("list", length(lff$fnfs))
  names(cmd) <- lff$fnfs
  output_R1 <- file.path(folder_output, basename(lff$fnfs))
  output_R2 <- if (paired_end && !single_file_output) {
    file.path(folder_output, basename(lff$fnrs))
  } else {
    NULL
  }

  for (i in seq_along(lff$fnfs)) {
    f_R1 <- lff$fnfs[i]
    base_name <- sub("\\.(fastq|fq)(\\.(gz|zst))?$", "", basename(f_R1))
    out_R1 <- output_R1[i]
    log_file <- file.path(folder_output, paste0(base_name, ".itsxpress.log"))

    itsxpress_cmd <- paste0(
      args_before_itsxpress,
      shQuote(itsxpresspath),
      " --fastq ",
      shQuote(f_R1),
      " --outfile ",
      shQuote(out_R1),
      " --region ",
      region,
      " --taxa ",
      taxa,
      " --cluster_id ",
      cluster_id,
      " --threads ",
      as.integer(nproc),
      " --log ",
      shQuote(log_file)
    )

    if (paired_end) {
      f_R2 <- lff$fnrs[i]
      itsxpress_cmd <- paste0(
        itsxpress_cmd,
        " --fastq2 ",
        shQuote(f_R2)
      )
      if (!single_file_output) {
        itsxpress_cmd <- paste0(
          itsxpress_cmd,
          " --outfile2 ",
          shQuote(output_R2[i])
        )
      }
    } else {
      itsxpress_cmd <- paste0(itsxpress_cmd, " --single_end")
    }

    if (reversed_primers) {
      itsxpress_cmd <- paste0(itsxpress_cmd, " --reversed_primers")
    }
    if (trim_ccs) {
      itsxpress_cmd <- paste0(itsxpress_cmd, " --trim-ccs")
    }
    if (nzchar(extra_itsxpress_args)) {
      itsxpress_cmd <- paste0(itsxpress_cmd, " ", extra_itsxpress_args)
    }

    cmd[[f_R1]] <- itsxpress_cmd
  }

  if (!cmd_is_run) {
    return(cmd)
  }

  if (
    !is_itsxpress_installed(
      path = itsxpresspath,
      args_before_itsxpress = args_before_itsxpress
    )
  ) {
    cli::cli_abort(c(
      "{.code itsxpress} was not found or could not be run.",
      "i" = "Install it with {.fn install_itsxpress} and pass the returned string to {.arg args_before_itsxpress}, or check {.fn is_itsxpress_installed}."
    ))
  }

  if (!dir.exists(folder_output)) {
    dir.create(folder_output, recursive = TRUE)
  }

  script_path <- tempfile("script_itsxpress_", fileext = ".sh")
  on.exit(unlink(script_path), add = TRUE)
  writeLines(unlist(cmd), script_path)
  status <- system2("bash", shQuote(script_path))
  if (status != 0) {
    cli::cli_abort(c(
      "ITSxpress failed (exit status {status}).",
      "i" = "Check the log files in {.path {folder_output}}."
    ))
  }
  message(
    "Output files are available in the folder ",
    normalizePath(folder_output)
  )

  if (return_file_path) {
    if (paired_end && !single_file_output) {
      return(list("fnfs" = output_R1, "fnrs" = output_R2))
    } else {
      return(list("fnfs" = output_R1))
    }
  }

  invisible(NULL)
}
