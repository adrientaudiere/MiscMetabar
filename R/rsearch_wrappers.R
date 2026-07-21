################################################################################
# Thin wrappers around Rsearch functions covering capabilities MiscMetabar
# does not implement natively: read-quality visualisation from fastq files,
# standalone dereplication and reference-based chimera detection. De-novo
# chimera detection is intentionally NOT wrapped: see chimera_detection_vs().
################################################################################

#' Abort when the Rsearch package is not available
#' @param fn Name of the calling function, used in the error message.
#' @keywords internal
#' @noRd
check_rsearch <- function(fn) {
  if (!requireNamespace("Rsearch", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "{.fn {fn}} requires the {.pkg Rsearch} package.",
        i = 'Install it with {.code install.packages("Rsearch")}.'
      )
    )
  }
}

################################################################################
#' Plot the quality of reads in a fastq file
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Wrapper of [Rsearch::plot_read_quality()]. Plot the per-base quality score
#' distribution of the reads of a fastq file, either as raw Phred scores or
#' converted to expected-error (EE) rates.
#'
#' @param fastq_input (required) Path to a fastq file (optionally gzipped).
#' @param use_ee_rate (logical, default FALSE) Convert quality scores to
#'   expected error rates.
#' @param plot_title (logical, default TRUE) Add a title to the plot.
#' @param alpha (numeric, default 0.5) Transparency of the quality lines.
#'
#' @returns A ggplot2-based plot object.
#' @export
#' @author Adrien Taudière
#' @seealso [Rsearch::plot_read_quality()], [plot_ee_rate_dist()]
#'
#' @examplesIf requireNamespace("Rsearch", quietly = TRUE)
#' plot_read_quality(
#'   system.file("extdata", "ex.fastq", package = "MiscMetabar")
#' )
plot_read_quality <- function(
  fastq_input,
  use_ee_rate = FALSE,
  plot_title = TRUE,
  alpha = 0.5
) {
  check_rsearch("plot_read_quality")
  Rsearch::plot_read_quality(
    fastq_input = fastq_input,
    use_ee_rate = use_ee_rate,
    plot_title = plot_title,
    alpha = alpha
  )
}
################################################################################

################################################################################
#' Plot the distribution of expected error (EE) rates of reads
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Wrapper of [Rsearch::plot_ee_rate_dist()]. Plot the distribution of the
#' expected error (EE) rate of the reads of a fastq file, a key diagnostic
#' before quality filtering (e.g. with `maxEE` thresholds).
#'
#' @param fastq_input (required) Path to a fastq file (optionally gzipped).
#' @param n_bins (integer, default 30) Number of histogram bins.
#' @param plot_title (character) Title of the plot.
#'
#' @returns A [ggplot2::ggplot()] object.
#' @export
#' @author Adrien Taudière
#' @seealso [Rsearch::plot_ee_rate_dist()], [plot_read_quality()]
#'
#' @examplesIf requireNamespace("Rsearch", quietly = TRUE)
#' plot_ee_rate_dist(
#'   system.file("extdata", "ex.fastq", package = "MiscMetabar")
#' )
plot_ee_rate_dist <- function(
  fastq_input,
  n_bins = 30,
  plot_title = "Distribution of the expected error (EE) rate of reads"
) {
  check_rsearch("plot_ee_rate_dist")
  Rsearch::plot_ee_rate_dist(
    fastq_input = fastq_input,
    n_bins = n_bins,
    plot_title = plot_title
  )
}
################################################################################

################################################################################
#' Dereplicate sequences in a fastx file
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Wrapper of [Rsearch::vs_fastx_uniques()] (vsearch `--fastx_uniques`).
#' Dereplicate the reads of a fasta or fastq file: identical sequences are
#' merged and their abundances recorded as `;size=N` annotations in the
#' output headers.
#'
#' @param fastx_input (required) Path to a fasta or fastq file (optionally
#'   gzipped).
#' @param output_format (default `"fasta"`) Output format: `"fasta"` or
#'   `"fastq"`.
#' @param fastx_output Path to the output file. A temporary file is used
#'   when `NULL`.
#' @param minuniquesize (integer, default 1) Discard unique sequences with
#'   an abundance lower than this value.
#' @param strand (default `"plus"`) Compare sequences on the `"plus"` strand
#'   only or on `"both"` strands.
#' @param sizein (logical, default TRUE) Take `;size=N` annotations of the
#'   input into account.
#' @param sizeout (logical, default TRUE) Add `;size=N` annotations to the
#'   output.
#' @param relabel Relabel output sequences with this prefix followed by a
#'   sequential number.
#' @param relabel_sha1 (logical, default FALSE) Relabel output sequences
#'   with their SHA1 digest.
#' @param fastq_qout_max (logical, default FALSE) Report the maximum quality
#'   score in fastq output.
#' @param fasta_width (integer, default 0) Width of fasta sequence lines
#'   (0 for no wrapping).
#' @param sample Optional sample name annotation added to output headers.
#' @param vsearch_options Additional vsearch command-line options.
#' @param tmpdir Directory for temporary files.
#'
#' @returns `NULL` (invisibly). The dereplicated sequences are written to
#'   `fastx_output`.
#' @export
#' @author Adrien Taudière
#' @seealso [Rsearch::vs_fastx_uniques()], [vs_uchime_ref()], [reads2otu()]
#'
#' @examplesIf requireNamespace("Rsearch", quietly = TRUE) && MiscMetabar::is_vsearch_installed()
#' vs_fastx_uniques(
#'   system.file("extdata", "ex_little.fasta", package = "MiscMetabar"),
#'   fastx_output = tempfile(fileext = ".fasta")
#' )
vs_fastx_uniques <- function(
  fastx_input,
  output_format = "fasta",
  fastx_output = NULL,
  minuniquesize = 1,
  strand = "plus",
  sizein = TRUE,
  sizeout = TRUE,
  relabel = NULL,
  relabel_sha1 = FALSE,
  fastq_qout_max = FALSE,
  fasta_width = 0,
  sample = NULL,
  vsearch_options = NULL,
  tmpdir = NULL
) {
  check_rsearch("vs_fastx_uniques")
  Rsearch::vs_fastx_uniques(
    fastx_input = fastx_input,
    output_format = output_format,
    fastx_output = fastx_output,
    minuniquesize = minuniquesize,
    strand = strand,
    sizein = sizein,
    sizeout = sizeout,
    relabel = relabel,
    relabel_sha1 = relabel_sha1,
    fastq_qout_max = fastq_qout_max,
    fasta_width = fasta_width,
    sample = sample,
    vsearch_options = vsearch_options,
    tmpdir = tmpdir
  )
}
################################################################################

################################################################################
#' Reference-based chimera detection with vsearch
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Wrapper of [Rsearch::vs_uchime_ref()] (vsearch `--uchime_ref`). Detect
#' chimeric sequences in a fasta file by comparison against a reference
#' database. Complements [chimera_detection_vs()], which performs *de novo*
#' detection without a reference database.
#'
#' @param fasta_input (required) Path to the query fasta file.
#' @param database (required) Path to the reference database fasta file.
#' @param nonchimeras Path to the output fasta file of non-chimeric
#'   sequences. A temporary file is used when `NULL`.
#' @param chimeras Path to the output fasta file of chimeric sequences.
#'   A temporary file is used when `NULL`.
#' @param sizein (logical, default TRUE) Take `;size=N` annotations of the
#'   input into account.
#' @param sizeout (logical, default TRUE) Add `;size=N` annotations to the
#'   output.
#' @param relabel Relabel output sequences with this prefix followed by a
#'   sequential number.
#' @param relabel_sha1 (logical, default FALSE) Relabel output sequences
#'   with their SHA1 digest.
#' @param fasta_width (integer, default 0) Width of fasta sequence lines
#'   (0 for no wrapping).
#' @param sample Optional sample name annotation added to output headers.
#' @param log_file Path to a log file.
#' @param threads (integer, default 1) Number of threads.
#' @param vsearch_options Additional vsearch command-line options.
#' @param tmpdir Directory for temporary files.
#'
#' @returns `NULL` (invisibly). The non-chimeric and chimeric sequences
#'   are written to the `nonchimeras` and `chimeras` files, respectively.
#' @export
#' @author Adrien Taudière
#' @seealso [Rsearch::vs_uchime_ref()], [chimera_detection_vs()],
#'   [chimera_removal_vs()]
#'
#' @examplesIf requireNamespace("Rsearch", quietly = TRUE) && MiscMetabar::is_vsearch_installed()
#' vs_uchime_ref(
#'   system.file("extdata", "ex.fasta", package = "MiscMetabar"),
#'   database = system.file(
#'     "extdata",
#'     "100_sp_UNITE_sh_general_release_dynamic.fasta",
#'     package = "MiscMetabar"
#'   )
#' )
vs_uchime_ref <- function(
  fasta_input,
  database,
  nonchimeras = NULL,
  chimeras = NULL,
  sizein = TRUE,
  sizeout = TRUE,
  relabel = NULL,
  relabel_sha1 = FALSE,
  fasta_width = 0,
  sample = NULL,
  log_file = NULL,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
) {
  check_rsearch("vs_uchime_ref")
  Rsearch::vs_uchime_ref(
    fasta_input = fasta_input,
    database = database,
    nonchimeras = nonchimeras,
    chimeras = chimeras,
    sizein = sizein,
    sizeout = sizeout,
    relabel = relabel,
    relabel_sha1 = relabel_sha1,
    fasta_width = fasta_width,
    sample = sample,
    log_file = log_file,
    threads = threads,
    vsearch_options = vsearch_options,
    tmpdir = tmpdir
  )
}
################################################################################
