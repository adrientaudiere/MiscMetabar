#' Align the reference sequences of a phyloseq object
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Build a multiple sequence alignment from the `refseq` slot of a
#' \code{\link[phyloseq]{phyloseq-class}} object, or from a set of sequences
#' given directly. Two backends are available:
#'
#' * `"decipher"` (default), the pure-R [DECIPHER::AlignSeqs()], which needs no
#'   external program;
#' * `"mafft"`, the MAFFT command-line aligner driven through [ips::mafft()],
#'   which is considerably faster on large sets of sequences and is the usual
#'   choice for metabarcoding-scale data.
#'
#' The alignment is returned as a `DNAStringSet`, **not** as a phyloseq object.
#' Gap characters are valid DNA letters, so an alignment can technically be
#' written back into a `refseq` slot, but the sequences of that slot are
#' consumed as unaligned elsewhere in the pqverse (clustering, BLAST, primer
#' trimming), so MiscMetabar never does it for you. Feed the returned alignment
#' to [build_phytree_pq()] or [Biostrings::writeXStringSet()] instead.
#'
#' @param x (required) A \code{\link[phyloseq]{phyloseq-class}} object with a
#'   non-empty `refseq` slot, a `DNAStringSet` (or any `XStringSet`), or a
#'   `DNAbin` object.
#' @param method One of `"decipher"` (default) or `"mafft"`.
#' @param exec Path to the MAFFT executable. Only used when
#'   `method = "mafft"`. Default to NULL, in which case MAFFT is looked up with
#'   [is_mafft_installed()]: first the `MiscMetabar.mafftpath` option, then the
#'   system `PATH`.
#' @param mafft_method The MAFFT strategy passed to [ips::mafft()], e.g.
#'   `"auto"` (default, lets MAFFT pick from the size of the input),
#'   `"localpair"` (L-INS-i), `"globalpair"` (G-INS-i) or `"retree 2"`.
#' @param thread Number of threads given to MAFFT. Default to -1, i.e. MAFFT
#'   detects the number of available cores itself.
#' @param force Logical, if TRUE sequences that are already all of the same
#'   length are realigned anyway. Default to FALSE, in which case such
#'   sequences are assumed to be aligned already and returned untouched.
#' @param verbose Logical, if TRUE report the aligner used and the width of the
#'   resulting alignment. Default to FALSE.
#' @param ... Additional arguments passed on to [DECIPHER::AlignSeqs()] or
#'   [ips::mafft()], depending on `method`.
#'
#' @return A `DNAStringSet` of aligned sequences, all of the same width, in the
#'   order of the input.
#'
#' @section Installing MAFFT:
#' MAFFT is an external program and is not installed by MiscMetabar. It is
#' available from \url{https://mafft.cbrc.jp/alignment/software/}, and from the
#' usual package managers:
#'
#' * Debian / Ubuntu: `sudo apt install mafft`
#' * macOS (Homebrew): `brew install mafft`
#' * conda: `conda install -c bioconda mafft`
#' * Windows: use the installer from the MAFFT website, or run it under WSL.
#'
#' If the executable is not on the `PATH`, either pass its path with `exec` or
#' set it once for the session with
#' `options(MiscMetabar.mafftpath = "/path/to/mafft")`.
#'
#' @author Adrien Taudière
#'
#' @seealso [is_mafft_installed()], [build_phytree_pq()]
#'
#' @examples
#' \donttest{
#' data(data_fungi_mini)
#' df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 12000)
#'
#' # Pure-R alignment, no external program required
#' ali <- align_pq(df, method = "decipher", verbose = TRUE)
#' ali
#' }
#'
#' \dontrun{
#' # The same alignment with MAFFT, much faster on large refseq slots
#' ali_mafft <- align_pq(df, method = "mafft", verbose = TRUE)
#'
#' # A more accurate (and slower) MAFFT strategy
#' ali_linsi <- align_pq(df, method = "mafft", mafft_method = "localpair")
#'
#' # An explicit path to the executable
#' ali <- align_pq(df, method = "mafft", exec = "/usr/local/bin/mafft")
#'
#' # Sequences can also be given directly
#' align_pq(phyloseq::refseq(df), method = "mafft")
#'
#' # Feed the alignment to a tree builder or write it to a FASTA file
#' Biostrings::writeXStringSet(ali_mafft, "refseq_aligned.fasta")
#' }
#' @export
align_pq <- function(
  x,
  method = c("decipher", "mafft"),
  exec = NULL,
  mafft_method = "auto",
  thread = -1,
  force = FALSE,
  verbose = FALSE,
  ...
) {
  method <- match.arg(method)
  dna <- as_dna_stringset(x)

  if (length(dna) < 2) {
    cli::cli_abort(
      "At least two sequences are required to build an alignment, not
       {length(dna)}."
    )
  }

  if (!force && length(unique(Biostrings::width(dna))) == 1) {
    if (verbose) {
      cli::cli_inform(c(
        "i" = "The {length(dna)} sequences already share a width of
               {unique(Biostrings::width(dna))}; returned untouched.",
        "i" = "Use {.code force = TRUE} to realign them."
      ))
    }
    return(dna)
  }

  ali <- if (method == "mafft") {
    align_with_mafft(dna, exec, mafft_method, thread, verbose, ...)
  } else {
    align_with_decipher(dna, verbose, ...)
  }

  if (verbose) {
    cli::cli_inform(c(
      "v" = "{length(ali)} sequences aligned with {.field {method}}
             over {unique(Biostrings::width(ali))} positions."
    ))
  }

  return(ali)
}

#' Is MAFFT available?
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Check whether the MAFFT command-line aligner used by
#' `align_pq(method = "mafft")` is installed. Useful to guard examples, tests
#' and vignette chunks.
#'
#' @param path Optional path to the MAFFT executable. Default to NULL, in which
#'   case MAFFT is looked up in two places, in order: the
#'   `MiscMetabar.mafftpath` option, then the system `PATH`.
#'
#' @return A logical of length one. FALSE when the `ips` package is not
#'   installed, so that the check also guards the R-level dependency.
#'
#' @author Adrien Taudière
#'
#' @seealso [align_pq()], [is_vsearch_installed()]
#'
#' @examples
#' is_mafft_installed()
#' @export
is_mafft_installed <- function(path = NULL) {
  if (!requireNamespace("ips", quietly = TRUE)) {
    return(FALSE)
  }
  if (!is.null(path)) {
    return(file.exists(path))
  }
  exec <- find_mafft()
  nzchar(exec) && file.exists(exec)
}

#' Locate the MAFFT executable
#'
#' Resolution order: the `MiscMetabar.mafftpath` option, then the system
#' `PATH`. MAFFT is a system package rather than a program MiscMetabar
#' compiles, so there is no user-data-directory step as in [find_vsearch()].
#'
#' @return A length-one character path, or `""` when nothing was found.
#' @noRd
#' @keywords internal
find_mafft <- function() {
  opt <- getOption("MiscMetabar.mafftpath")
  if (!is.null(opt) && nzchar(opt)) {
    return(opt)
  }
  unname(Sys.which("mafft"))
}

#' Coerce the `x` argument of [align_pq()] to a `DNAStringSet`
#'
#' @inheritParams align_pq
#' @return A `DNAStringSet`.
#' @noRd
#' @keywords internal
as_dna_stringset <- function(x) {
  if (inherits(x, "phyloseq")) {
    if (is.null(x@refseq)) {
      cli::cli_abort(c(
        "The {.arg x} object has no {.field refseq} slot.",
        "i" = "{.fn align_pq} aligns reference sequences."
      ))
    }
    return(Biostrings::DNAStringSet(x@refseq))
  }

  if (inherits(x, "DNAbin")) {
    return(dnabin_to_stringset(x))
  }

  if (methods::is(x, "XStringSet")) {
    return(Biostrings::DNAStringSet(x))
  }

  cli::cli_abort(
    "{.arg x} must be a {.cls phyloseq}, {.cls XStringSet} or {.cls DNAbin}
     object, not {.cls {class(x)}}."
  )
}

#' Align with [DECIPHER::AlignSeqs()]
#'
#' @inheritParams align_pq
#' @param dna A `DNAStringSet`.
#' @return An aligned `DNAStringSet`.
#' @noRd
#' @keywords internal
align_with_decipher <- function(dna, verbose = FALSE, ...) {
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg DECIPHER} is required for
       {.code method = \"decipher\"}.",
      "i" = "Install it with
             {.code BiocManager::install(\"DECIPHER\")}, or use
             {.code method = \"mafft\"}."
    ))
  }
  DECIPHER::AlignSeqs(dna, anchor = NA, verbose = verbose, ...)
}

#' Align with MAFFT through [ips::mafft()]
#'
#' `ips::mafft()` works on `DNAbin` objects and returns a `DNAbin` matrix, so
#' the sequences make a round trip. Sequence names survive it, but the letters
#' come back lower-case, hence the `toupper()`.
#'
#' @inheritParams align_pq
#' @param dna A `DNAStringSet`.
#' @return An aligned `DNAStringSet`.
#' @noRd
#' @keywords internal
align_with_mafft <- function(
  dna,
  exec = NULL,
  mafft_method = "auto",
  thread = -1,
  verbose = FALSE,
  ...
) {
  if (!requireNamespace("ips", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg ips} is required for {.code method = \"mafft\"}.",
      "i" = "Install it with {.code install.packages(\"ips\")}."
    ))
  }
  exec <- resolve_mafft_exec(exec)

  ali <- ips::mafft(
    ape::as.DNAbin(dna),
    method = mafft_method,
    exec = exec,
    thread = thread,
    quiet = !verbose,
    ...
  )

  if (!inherits(ali, "DNAbin") || is.null(dim(ali))) {
    cli::cli_abort(c(
      "{.field mafft} returned no alignment.",
      "i" = "Check that {.path {exec}} runs from a terminal.",
      "i" = "Use {.code verbose = TRUE} to see its output."
    ))
  }

  dnabin_to_stringset(ali)
}

#' Resolve the path to the MAFFT executable
#'
#' @inheritParams align_pq
#' @return A length-one character path to the executable.
#' @noRd
#' @keywords internal
resolve_mafft_exec <- function(exec) {
  if (is.null(exec)) {
    exec <- find_mafft()
  }
  if (!nzchar(exec) || !file.exists(exec)) {
    cli::cli_abort(c(
      "The {.field mafft} executable was not found.",
      "i" = "Install it, e.g. {.code sudo apt install mafft},
             {.code brew install mafft} or
             {.code conda install -c bioconda mafft}.",
      "i" = "Give its path with {.arg exec}, or set it once with
             {.code options(MiscMetabar.mafftpath = \"/path/to/mafft\")}.",
      "i" = "Installation instructions:
             {.url https://mafft.cbrc.jp/alignment/software/}."
    ))
  }
  exec
}

#' Convert a `DNAbin` alignment back to a `DNAStringSet`
#'
#' @param x A `DNAbin` matrix or list.
#' @return A `DNAStringSet`.
#' @noRd
#' @keywords internal
dnabin_to_stringset <- function(x) {
  chars <- as.character(x)
  seqs <- if (is.matrix(chars)) {
    apply(chars, 1, paste, collapse = "")
  } else {
    vapply(chars, paste, character(1), collapse = "")
  }
  Biostrings::DNAStringSet(toupper(seqs))
}
