################################################################################
#' Find the ITSx executable
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Resolution order: the `MiscMetabar.itsxpath` option (if set), then the
#' system `PATH`. ITSx is usually installed in a conda environment (see
#' [install_itsx()]); in that case it is not on the `PATH` of R and the
#' environment is activated with the `args_before_itsx` argument of
#' [itsx_pq()] and [is_itsx_installed()] instead.
#'
#' @return A character string with the path to ITSx, or `"ITSx"` as a
#'   fallback (relying on `PATH` resolution, possibly after
#'   `args_before_itsx`).
#' @export
#' @examples
#' find_itsx()
#' @author Adrien Taudière
#' @seealso [itsx_pq()], [is_itsx_installed()], [install_itsx()]
find_itsx <- function() {
  opt <- getOption("MiscMetabar.itsxpath")
  if (!is.null(opt) && nzchar(opt)) {
    return(opt)
  }
  on_path <- unname(Sys.which("ITSx"))
  if (nzchar(on_path)) {
    return(on_path)
  }
  "ITSx"
}

################################################################################
#' Test if ITSx is installed
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage.
#'
#' @param path (default: [find_itsx()]) Path to ITSx.
#' @param args_before_itsx (String, default "") A one line bash command run
#'   before ITSx, e.g. the value returned by [install_itsx()]:
#'   `"source ~/miniforge3/etc/profile.d/conda.sh && conda activate itsxenv && "`.
#' @export
#' @return A logical that says if ITSx can be run.
#'
#' @examples
#' MiscMetabar::is_itsx_installed()
#' @author Adrien Taudière
#' @seealso [find_itsx()], [install_itsx()], [itsx_pq()]
is_itsx_installed <- function(path = find_itsx(), args_before_itsx = "") {
  status <- suppressWarnings(system2(
    "bash",
    c("-c", shQuote(paste0(args_before_itsx, shQuote(path), " --help"))),
    stdout = FALSE,
    stderr = FALSE
  ))
  identical(as.integer(status), 0L)
}

################################################################################
#' Install ITSx in a conda environment
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' ITSx (Bengtsson-Palme et al. 2013) is a Perl program that needs HMMER 3.
#' The simplest way to get both is the bioconda package, which this function
#' installs in a dedicated conda environment. It returns the string to pass
#' as `args_before_itsx` to [itsx_pq()] and [is_itsx_installed()], which
#' activates that environment before each call.
#'
#' @param env_name (default: "itsxenv") Name of the conda environment.
#' @param conda (default: NULL) Path to `mamba` or `conda`. If NULL,
#'   `mamba` then `conda` are looked up on the `PATH`.
#' @param force (default: FALSE) If TRUE, the environment is created even if
#'   one of that name already exists (conda then replaces it).
#'
#' @return The `args_before_itsx` string activating the environment
#'   (invisibly).
#'
#' @section Without conda:
#' Download ITSx from \url{https://microbiology.se/software/itsx/}, install
#' HMMER 3 (\url{http://hmmer.org/}; `sudo apt install hmmer` on Debian /
#' Ubuntu, `brew install hmmer` on macOS), put the `ITSx` script and its
#' `ITSx_db` folder on the `PATH`, or set
#' `options(MiscMetabar.itsxpath = "/path/to/ITSx")`.
#'
#' @export
#' @examples
#' \dontrun{
#' prelude <- install_itsx()
#' is_itsx_installed(args_before_itsx = prelude)
#' }
#' @author Adrien Taudière
#' @seealso [itsx_pq()], [is_itsx_installed()], [find_itsx()]
install_itsx <- function(env_name = "itsxenv", conda = NULL, force = FALSE) {
  if (is.null(conda)) {
    candidates <- unname(Sys.which(c("mamba", "conda")))
    conda <- candidates[nzchar(candidates)][1]
  }
  if (is.na(conda) || !nzchar(conda)) {
    cli::cli_abort(c(
      "Neither {.code mamba} nor {.code conda} was found on the PATH.",
      "i" = "Install Miniforge (\\url{https://github.com/conda-forge/miniforge})
      or follow the section {.emph Without conda} of {.fn install_itsx}."
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

  if (!force && is_itsx_installed(args_before_itsx = prelude)) {
    message(
      "ITSx is already installed in the conda environment '",
      env_name,
      "'. Use force = TRUE to reinstall."
    )
    return(invisible(prelude))
  }

  message("Creating the conda environment '", env_name, "' with ITSx...")
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
      "itsx"
    )
  )
  if (status != 0 || !is_itsx_installed(args_before_itsx = prelude)) {
    cli::cli_abort(
      "The installation of ITSx in the conda environment {.val {env_name}} failed."
    )
  }
  message("ITSx installed. Use args_before_itsx = \"", prelude, "\"")
  invisible(prelude)
}

# Organisms of the ITSx profile sets, by the one-letter code ITSx writes in
# its fasta headers (the table of the ITSx 1.1.3 script).
itsx_origin_codes <- c(
  A = "Alveolata",
  B = "Bryophyta",
  C = "Bacillariophyta",
  D = "Amoebozoa",
  E = "Euglenozoa",
  F = "Fungi",
  G = "Chlorophyta",
  H = "Rhodophyta",
  I = "Phaeophyceae",
  L = "Marchantiophyta",
  M = "Metazoa",
  N = "Microsporidia",
  O = "Oomycota",
  P = "Haptophyceae",
  Q = "Raphidophyceae",
  R = "Rhizaria",
  S = "Synurophyceae",
  T = "Tracheophyta",
  U = "Eustigmatophyceae",
  X = "Apusozoa",
  Y = "Parabasalia"
)

# Split ITSx fasta headers into the id and the one-letter origin code. ITSx
# writes "<id>|<code>|<region> Extracted ..." in the region files (ITS1.fasta,
# SSU.fasta...) but "<id>|<code> fungi ITS sequence (503 bp) ..." in
# full.fasta. The id is taken greedily, so an id that itself contains "|" is
# kept whole.
#' @noRd
#' @keywords internal
parse_itsx_headers <- function(headers) {
  pattern <- "^(.*)\\|([A-Z.])(\\|[^ |]+)?( .*)?$"
  ok <- grepl(pattern, headers)
  data.frame(
    id = ifelse(ok, sub(pattern, "\\1", headers), sub(" .*$", "", headers)),
    code = ifelse(ok, sub(pattern, "\\2", headers), NA_character_),
    stringsAsFactors = FALSE
  )
}

# Run ITSx once on `seqs` and return every region it wrote, named like `seqs`,
# and the origin code of each detected sequence.
#' @noRd
#' @keywords internal
run_itsx_seqs <- function(
  seqs,
  organism_groups = "all",
  nproc = 1,
  itsxpath = find_itsx(),
  args_before_itsx = ""
) {
  work_dir <- tempfile("itsx_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  # Short ids: the names of `seqs` may hold "|" or spaces, which ITSx and its
  # headers would mangle.
  ids <- paste0("s", seq_along(seqs))
  input <- file.path(work_dir, "input.fasta")
  prefix <- file.path(work_dir, "itsx")
  seqs_out <- seqs
  names(seqs_out) <- ids
  Biostrings::writeXStringSet(seqs_out, input, width = 20000)

  cmd <- paste0(
    args_before_itsx,
    shQuote(itsxpath),
    " -i ",
    shQuote(input),
    " -o ",
    shQuote(prefix),
    " -t ",
    shQuote(organism_groups),
    " --cpu ",
    as.integer(nproc),
    " --preserve F --save_regions all --graphical F --silent T"
  )
  status <- system2("bash", c("-c", shQuote(cmd)), stdout = FALSE)
  if (status != 0) {
    cli::cli_abort(c(
      "ITSx failed (exit status {status}).",
      "i" = "Check {.fn is_itsx_installed} with the same {.arg args_before_itsx}."
    ))
  }

  region_files <- c(
    SSU = "SSU",
    ITS1 = "ITS1",
    `5.8S` = "5_8S",
    ITS2 = "ITS2",
    LSU = "LSU",
    full = "full"
  )
  regions <- lapply(region_files, function(r) {
    f <- paste0(prefix, ".", r, ".fasta")
    if (!file.exists(f) || file.size(f) == 0) {
      return(Biostrings::DNAStringSet())
    }
    Biostrings::readDNAStringSet(f)
  })
  origin <- stats::setNames(rep(NA_character_, length(ids)), ids)
  for (r in names(regions)) {
    if (length(regions[[r]]) == 0) {
      next
    }
    parsed <- parse_itsx_headers(names(regions[[r]]))
    known <- parsed$id %in% ids
    if (!all(known)) {
      cli::cli_abort(
        "Could not read the ITSx header {.val {names(regions[[r]])[!known][1]}}."
      )
    }
    origin[parsed$id] <- ifelse(
      is.na(origin[parsed$id]),
      parsed$code,
      origin[parsed$id]
    )
    names(regions[[r]]) <- names(seqs)[match(parsed$id, ids)]
  }
  names(origin) <- names(seqs)
  list(regions = regions, origin = origin)
}

################################################################################
#' Extract the ITS region of the sequences of a phyloseq object with ITSx
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Run ITSx (Bengtsson-Palme et al. 2013) on the `refseq` slot of a
#' phyloseq object and return the object with its `refseq` replaced by the
#' extracted region. ITSx locates the conserved flanks (end of the SSU, 5.8S,
#' start of the LSU) with HMM profiles, so the region is cut where the
#' ribosomal genes actually are, whatever the primers.
#'
#' Removing the conserved flanks matters beyond tidiness: they are nearly
#' identical in every reference record, so an alignment-based assignment such
#' as [assign_blastn()] aligns a query that still carries them against
#' almost the whole database. On fungal ITS2 amplicons, extracting the ITS2
#' made blastn more than ten times faster.
#'
#' ITSx also reports, for each sequence it detects, the organism group whose
#' profiles fit it best (its *putative origin*). `add_origin` keeps it as a
#' column of the `tax_table`, and `remove_other_origin` drops the taxa of
#' another origin, a cheap filter against non-target sequences.
#'
#' @inheritParams clean_pq
#' @param region (default: "ITS1") The region to keep: `"ITS1"`, `"ITS2"`
#'   or `"full"` (ITS1 + 5.8S + ITS2, only for the sequences in which ITSx
#'   finds all three).
#' @param organism_groups (default: "all") The ITSx profile sets to search,
#'   passed to `-t`: `"all"`, or comma-separated one-letter codes such as
#'   `"F"` (fungi) or `"F,T"`. See Details for the codes.
#' @param keep_undetected (logical, default TRUE) If TRUE, the taxa in which
#'   ITSx does not find `region` keep their full sequence, so the set of taxa
#'   is unchanged. If FALSE, they are removed.
#' @param add_origin (logical, default TRUE) If TRUE, a column `origin_col`
#'   is added to the `tax_table` with the putative origin ITSx gives to each
#'   taxon (e.g. `"Fungi"`), NA when ITSx detects nothing in the sequence.
#' @param origin_col (default: "ITSx_origin") Name of that column.
#' @param remove_other_origin (logical, default FALSE) If TRUE, the taxa
#'   whose putative origin is not in `keep_origin` are removed, including the
#'   taxa in which ITSx detects nothing.
#' @param keep_origin (default: "F") One-letter codes of the origins kept when
#'   `remove_other_origin` is TRUE.
#' @param duplicated_seqs (default: "merge") What to do with taxa whose
#'   extracted sequences are identical (they differed only in the removed
#'   flanks): `"merge"` merges them with [merge_taxa_vec()] into the most
#'   abundant one (its taxonomy is kept), `"remove"` keeps the most abundant
#'   and removes the others (ties: the first one).
#' @param nproc (default: 1) Number of CPUs given to ITSx (`--cpu`).
#' @param itsxpath (default: [find_itsx()]) Path to ITSx.
#' @param args_before_itsx (String, default "") A one line bash command run
#'   before ITSx, e.g. the conda activation returned by [install_itsx()].
#' @param verbose (logical, default TRUE) If TRUE, report how many taxa were
#'   detected, removed or merged.
#'
#' @details
#' Origin codes of ITSx 1.1.3: A Alveolata, B Bryophyta, C Bacillariophyta,
#' D Amoebozoa, E Euglenozoa, F Fungi, G Chlorophyta, H Rhodophyta,
#' I Phaeophyceae, L Marchantiophyta, M Metazoa, N Microsporidia, O Oomycota,
#' P Haptophyceae, Q Raphidophyceae, R Rhizaria, S Synurophyceae,
#' T Tracheophyta, U Eustigmatophyceae, X Apusozoa, Y Parabasalia.
#'
#' The putative origin is the profile set with the best score; on short or
#' partial sequences (e.g. an ITS1 amplicon ending inside the 5.8S) it can
#' be wrong, so check it before using `remove_other_origin`.
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object whose
#'   `refseq` holds the extracted region.
#' @export
#' @references
#' Bengtsson-Palme J, Ryberg M, Hartmann M, et al. (2013). Improved software
#' detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of
#' fungi and other eukaryotes for analysis of environmental sequencing data.
#' *Methods in Ecology and Evolution* 4: 914–919.
#' \doi{10.1111/2041-210X.12073}
#' @examples
#' \dontrun{
#' prelude <- install_itsx()
#' d_its1 <- itsx_pq(data_fungi_mini, region = "ITS1", args_before_itsx = prelude)
#' table(d_its1@tax_table[, "ITSx_origin"], useNA = "ifany")
#' }
#' @author Adrien Taudière
#' @seealso [install_itsx()], [is_itsx_installed()], [cutadapt_remove_primers()]
itsx_pq <- function(
  physeq,
  region = c("ITS1", "ITS2", "full"),
  organism_groups = "all",
  keep_undetected = TRUE,
  add_origin = TRUE,
  origin_col = "ITSx_origin",
  remove_other_origin = FALSE,
  keep_origin = "F",
  duplicated_seqs = c("merge", "remove"),
  nproc = 1,
  itsxpath = find_itsx(),
  args_before_itsx = "",
  verbose = TRUE
) {
  region <- match.arg(region)
  duplicated_seqs <- match.arg(duplicated_seqs)
  verify_pq(physeq)
  if (is.null(physeq@refseq)) {
    cli::cli_abort("{.arg physeq} has no {.field refseq} slot.")
  }

  res <- run_itsx_seqs(
    physeq@refseq,
    organism_groups = organism_groups,
    nproc = nproc,
    itsxpath = itsxpath,
    args_before_itsx = args_before_itsx
  )
  extracted <- res$regions[[region]]
  origin_code <- res$origin[phyloseq::taxa_names(physeq)]
  detected <- phyloseq::taxa_names(physeq) %in% names(extracted)

  keep <- rep(TRUE, phyloseq::ntaxa(physeq))
  if (!keep_undetected) {
    keep <- keep & detected
  }
  if (remove_other_origin) {
    keep <- keep & !is.na(origin_code) & origin_code %in% keep_origin
  }
  if (verbose) {
    message(
      sum(detected),
      " of ",
      length(detected),
      " taxa with ",
      region,
      " detected by ITSx; ",
      sum(!keep),
      " removed (",
      if (!keep_undetected) "undetected" else "",
      if (!keep_undetected && remove_other_origin) " or " else "",
      if (remove_other_origin) "other origin" else "",
      if (keep_undetected && !remove_other_origin) "none asked" else "",
      ")."
    )
  }

  new_seqs <- as.character(physeq@refseq)
  new_seqs[detected] <- as.character(extracted[
    phyloseq::taxa_names(physeq)[detected]
  ])
  new_physeq <- physeq
  new_physeq@refseq <- Biostrings::DNAStringSet(new_seqs)

  if (add_origin) {
    origin <- unname(itsx_origin_codes[origin_code])
    if (is.null(new_physeq@tax_table)) {
      tax <- matrix(
        origin,
        ncol = 1,
        dimnames = list(phyloseq::taxa_names(physeq), origin_col)
      )
    } else {
      tax <- as(new_physeq@tax_table, "matrix")
      tax <- cbind(tax, stats::setNames(origin, NULL))
      colnames(tax)[ncol(tax)] <- origin_col
    }
    new_physeq@tax_table <- phyloseq::tax_table(tax)
  }

  if (any(!keep)) {
    new_physeq <- phyloseq::prune_taxa(
      phyloseq::taxa_names(new_physeq)[keep],
      new_physeq
    )
  }

  seqs_now <- as.character(new_physeq@refseq)
  if (anyDuplicated(seqs_now)) {
    n_dup <- sum(duplicated(seqs_now))
    if (duplicated_seqs == "merge") {
      new_physeq <- merge_taxa_vec(
        new_physeq,
        group = match(seqs_now, unique(seqs_now)),
        tax_adjust = 0
      )
    } else {
      abundance <- phyloseq::taxa_sums(new_physeq)
      by_abundance <- order(-abundance, seq_along(seqs_now))
      dropped <- phyloseq::taxa_names(new_physeq)[by_abundance][
        duplicated(seqs_now[by_abundance])
      ]
      new_physeq <- phyloseq::prune_taxa(
        setdiff(phyloseq::taxa_names(new_physeq), dropped),
        new_physeq
      )
    }
    if (verbose) {
      message(
        n_dup,
        " taxa sharing their ",
        region,
        " with a more abundant taxon were ",
        if (duplicated_seqs == "merge") "merged into it." else "removed."
      )
    }
  }
  new_physeq
}
