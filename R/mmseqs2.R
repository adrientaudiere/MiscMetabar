################################################################################
#' Find the MMseqs2 binary
#'
#' @description
#' Looks for the MMseqs2 binary in three places, in order:
#' 1. The option `MiscMetabar.mmseqs2path` (if set).
#' 2. A local copy installed by [install_mmseqs2()] in the user data directory.
#' 3. The system `PATH`.
#'
#' @returns A character string with the path to the `mmseqs` binary.
#' @export
#' @author Adrien Taudière
#' @seealso [install_mmseqs2()], [is_mmseqs2_installed()], [assign_mmseqs2()]
find_mmseqs2 <- function() {
  # 1. User-set option
  opt <- getOption("MiscMetabar.mmseqs2path")
  if (!is.null(opt) && nzchar(opt)) {
    return(opt)
  }

  # 2. Installed copy in user data dir
  data_dir <- tools::R_user_dir("MiscMetabar", "data")
  bin_name <- if (.Platform$OS.type == "windows") "mmseqs.exe" else "mmseqs"
  local_bin <- file.path(data_dir, "bin", bin_name)
  if (file.exists(local_bin)) {
    return(local_bin)
  }

  # 3. Fallback to PATH
  "mmseqs"
}

################################################################################
#' Check whether MMseqs2 is installed and callable
#'
#' @description
#' Tries to run `mmseqs version` and returns `TRUE` if it succeeds.
#'
#' @param path Path to the `mmseqs` binary (default: [find_mmseqs2()]).
#' @returns Logical.
#' @export
#' @author Adrien Taudière
#' @seealso [find_mmseqs2()], [install_mmseqs2()], [assign_mmseqs2()]
#' @examples
#' is_mmseqs2_installed()
is_mmseqs2_installed <- function(path = find_mmseqs2()) {
  return(
    !inherits(
      try(
        system2(path, "version", stdout = TRUE, stderr = TRUE),
        silent = TRUE
      ),
      "try-error"
    )
  )
}

################################################################################
#' Install MMseqs2 from GitHub releases
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Downloads a pre-compiled MMseqs2 binary from
#' <https://mmseqs.com/latest/> and places it in the user data
#' directory for this package. Subsequent calls to [find_mmseqs2()] will
#' find it automatically.
#'
#' @param version Character. Either `"latest"` (default) or a specific
#'   release tag (e.g. `"17-b804f"`).
#' @param path Destination directory (default:
#'   `tools::R_user_dir("MiscMetabar", "data")`).
#' @param force Logical. Re-download even if already installed?
#' @returns The path to the installed binary (invisibly).
#' @export
#' @author Adrien Taudière
#' @seealso [find_mmseqs2()], [is_mmseqs2_installed()], [assign_mmseqs2()]
#' @examples
#' \dontrun{
#' install_mmseqs2()
#' install_mmseqs2(force = TRUE)
#' }
install_mmseqs2 <- function(
  version = "latest",
  path = tools::R_user_dir("MiscMetabar", "data"),
  force = FALSE
) {
  bin_name <- if (.Platform$OS.type == "windows") "mmseqs.exe" else "mmseqs"
  dest_bin <- file.path(path, "bin", bin_name)

  if (file.exists(dest_bin) && !force) {
    message("MMseqs2 is already installed at ", dest_bin)
    message("Use force = TRUE to reinstall.")
    return(invisible(dest_bin))
  }

  os <- tolower(Sys.info()[["sysname"]])

  # Determine platform tag for download URL
  platform_info <- switch(
    os,
    "darwin" = list(tag = "osx-universal", ext = ".tar.gz"),
    "linux" = {
      # Pick the best SIMD variant available
      cpu_flags <- tryCatch(
        readLines(file.path("", "proc", "cpuinfo")),
        error = function(e) character(0)
      )
      cpu_text <- paste(cpu_flags, collapse = " ")
      variant <- if (grepl("avx2", cpu_text, ignore.case = TRUE)) {
        "linux-avx2"
      } else if (grepl("sse4_1", cpu_text, ignore.case = TRUE)) {
        "linux-sse41"
      } else {
        "linux-sse2"
      }
      list(tag = variant, ext = ".tar.gz")
    },
    "windows" = list(tag = "win64", ext = ".zip"),
    stop("Unsupported operating system: ", os)
  )
  platform_tag <- platform_info$tag
  ext <- platform_info$ext

  archive_name <- paste0("mmseqs-", platform_tag, ext)
  download_url <- paste0("https://mmseqs.com/latest/", archive_name)

  message("Downloading MMseqs2 for ", platform_tag, " ...")
  archive_file <- file.path(tempdir(), archive_name)
  utils::download.file(download_url, archive_file, mode = "wb", quiet = TRUE)

  # Extract
  extract_dir <- file.path(tempdir(), "mmseqs2_extract")
  if (dir.exists(extract_dir)) {
    unlink(extract_dir, recursive = TRUE)
  }
  dir.create(extract_dir)

  if (ext == ".zip") {
    utils::unzip(archive_file, exdir = extract_dir)
  } else {
    utils::untar(archive_file, exdir = extract_dir)
  }

  # Find the binary inside the extracted archive
  extracted_bin <- list.files(
    extract_dir,
    pattern = paste0("^", bin_name, "$"),
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(extracted_bin) == 0) {
    stop("Could not find mmseqs binary in the downloaded archive.")
  }

  # Copy binary to destination

  dest_bin_dir <- file.path(path, "bin")
  dir.create(dest_bin_dir, recursive = TRUE, showWarnings = FALSE)

  file.copy(extracted_bin[[1]], dest_bin, overwrite = TRUE)
  if (.Platform$OS.type != "windows") {
    Sys.chmod(dest_bin, mode = "0755")
  }

  # On Windows, also copy DLLs
  if (os == "windows") {
    dlls <- list.files(
      dirname(extracted_bin[[1]]),
      pattern = "\\.dll$",
      full.names = TRUE
    )
    for (dll in dlls) {
      file.copy(dll, file.path(dest_bin_dir, basename(dll)), overwrite = TRUE)
    }
  }

  # Cleanup
  unlink(archive_file)
  unlink(extract_dir, recursive = TRUE)

  message("MMseqs2 installed at ", dest_bin)
  invisible(dest_bin)
}

################################################################################
# Build a synthetic MMseqs2 taxonomy database from SINTAX-formatted headers.
# Creates NCBI-format taxdump files and a sequence+taxonomy database so that
# easy-taxonomy (with --lca-mode) can be used on a plain SINTAX FASTA.
# Returns the path to the resulting MMseqs2 sequence database.
# @noRd
.build_sintax_mmseqs2_taxdb <- function(
  ref_ids,
  tax_strings,
  prefix_to_rank,
  ref_fasta_path,
  mmseqs2path,
  verbose = FALSE
) {
  # Build taxonomy tree -------------------------------------------------------
  nodes_env <- new.env(hash = TRUE, parent = emptyenv())
  next_id <- 2L
  assign(
    "root",
    list(taxid = 1L, parent = 1L, rank = "no rank", name = "root"),
    envir = nodes_env
  )
  get_or_create <- function(key, name, rank, parent_taxid) {
    node <- tryCatch(get(key, envir = nodes_env), error = \(e) NULL)
    if (is.null(node)) {
      id <- next_id
      next_id <<- next_id + 1L
      assign(
        key,
        list(taxid = id, parent = parent_taxid, rank = rank, name = name),
        envir = nodes_env
      )
      id
    } else {
      node$taxid
    }
  }

  seq_taxids <- integer(length(ref_ids))
  for (i in seq_along(ref_ids)) {
    parent_id <- 1L
    if (nzchar(tax_strings[i])) {
      for (part in strsplit(tax_strings[i], ",")[[1]]) {
        prefix <- substr(part, 1, 1)
        name <- sub("^[a-z]:", "", trimws(part))
        rank_name <- prefix_to_rank[prefix]
        if (!is.na(rank_name) && nzchar(name)) {
          key <- paste0(rank_name, ":", name)
          parent_id <- get_or_create(key, name, rank_name, parent_id)
        }
      }
    }
    leaf_id <- next_id
    next_id <- next_id + 1L
    assign(
      paste0("leaf:", ref_ids[i]),
      list(
        taxid = leaf_id,
        parent = parent_id,
        rank = "no rank",
        name = ref_ids[i]
      ),
      envir = nodes_env
    )
    seq_taxids[i] <- leaf_id
  }

  # Collect nodes and sort by taxid -------------------------------------------
  all_nodes <- as.list(nodes_env)
  tids <- vapply(all_nodes, \(n) n$taxid, integer(1))
  pids <- vapply(all_nodes, \(n) n$parent, integer(1))
  rnks <- vapply(all_nodes, \(n) n$rank, character(1))
  nms <- vapply(all_nodes, \(n) n$name, character(1))
  ord <- order(tids)
  tids <- tids[ord]
  pids <- pids[ord]
  rnks <- rnks[ord]
  nms <- nms[ord]

  # Write NCBI taxdump --------------------------------------------------------
  taxdump_dir <- file.path(tempdir(), "mmseqs2_taxdump")
  dir.create(taxdump_dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(
    paste0(tids, "\t|\t", nms, "\t|\t\t|\tscientific name\t|"),
    file.path(taxdump_dir, "names.dmp")
  )
  writeLines(
    paste0(
      tids,
      "\t|\t",
      pids,
      "\t|\t",
      rnks,
      "\t|\t\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t"
    ),
    file.path(taxdump_dir, "nodes.dmp")
  )
  writeLines(character(0), file.path(taxdump_dir, "merged.dmp"))

  # Write seqid -> taxid mapping ----------------------------------------------
  mapping_file <- file.path(tempdir(), "mmseqs2_seqtaxid.tsv")
  writeLines(paste0(ref_ids, "\t", seq_taxids), mapping_file)

  # Build MMseqs2 sequence database ------------------------------------------
  ref_db <- file.path(tempdir(), "mmseqs2_refDB")
  system2(
    mmseqs2path,
    paste("createdb", ref_fasta_path, ref_db),
    stdout = if (verbose) "" else TRUE,
    stderr = if (verbose) "" else TRUE
  )

  # Add taxonomy to the database ---------------------------------------------
  createtaxdb_tmp <- file.path(tempdir(), "mmseqs2_createtaxdb_tmp")
  dir.create(createtaxdb_tmp, showWarnings = FALSE, recursive = TRUE)
  system2(
    mmseqs2path,
    paste(
      "createtaxdb",
      ref_db,
      createtaxdb_tmp,
      "--ncbi-tax-dump",
      taxdump_dir,
      "--tax-mapping-file",
      mapping_file
    ),
    stdout = if (verbose) "" else TRUE,
    stderr = if (verbose) "" else TRUE
  )

  ref_db
}

################################################################################
#' Assign taxonomy using MMseqs2
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Use the [MMseqs2](https://github.com/soedinglab/MMseqs2) software to
#' assign taxonomy to sequences.
#'
#' The preferred usage is to provide a reference FASTA file in SINTAX format
#' via `ref_fasta`. The function builds a temporary MMseqs2 taxonomy database
#' from the SINTAX headers and then runs `mmseqs easy-taxonomy` with the
#' requested `--lca-mode`, giving the same LCA behaviour as the `database`
#' path.
#'
#' Alternatively, a pre-built MMseqs2 database with NCBI taxonomy can be
#' passed via the `database` parameter (created via `mmseqs createdb` +
#' `mmseqs createtaxdb`, or downloaded with `mmseqs databases`). In this
#' case, the MMseqs2 native `easy-taxonomy` LCA workflow is used. See the
#' [MMseqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki) for details.
#'
#' @inheritParams clean_pq
#' @param ref_fasta Either a [Biostrings::DNAStringSet] object or a path
#'   to a FASTA file in SINTAX format (taxonomy in headers after
#'   `;tax=`). Only used if `database` is not set.
#'   See [assign_sintax()] for the SINTAX format specification.
#' @param seq2search (optional) A [Biostrings::DNAStringSet] object. Use
#'   instead of `physeq` to search arbitrary sequences. Cannot be used
#'   together with `physeq`.
#' @param database (optional) Path to a pre-built MMseqs2 database with
#'   NCBI taxonomy information. Only used if `ref_fasta` is not set.
#' @param mmseqs2path Path to the `mmseqs` binary
#'   (default: [find_mmseqs2()]).
#' @param behavior Either `"return_matrix"` (default) or
#'   `"add_to_phyloseq"`:
#'   - `"return_matrix"`: return a data frame with taxonomic assignments.
#'   - `"add_to_phyloseq"`: return a phyloseq object with the
#'     taxonomy appended to the `tax_table` slot.
#' @param suffix (character) Suffix appended to new taxonomy column names
#'   (default: `"_mmseqs2"`).
#' @param lca_mode (integer) The LCA mode used by MMseqs2:
#'   - `1`: single search LCA
#'   - `3` (default): approximate 2bLCA (fast, recommended)
#'   - `4`: top-hit LCA (all equal-scoring top hits)
#' @param lca_ranks Character vector of NCBI taxonomy rank names passed
#'   to `--lca-ranks` (default:
#'   `c("superkingdom", "phylum", "class", "order", "family", "genus",
#'   "species")`).
#' @param column_names Character vector of output column names, must be
#'   the same length as `lca_ranks` (default:
#'   `c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
#'   "Species")`).
#' @param search_type (integer) MMseqs2 search type:
#'   - `0`: auto-detect
#'   - `2`: translated nucleotide
#'   - `3` (default): nucleotide
#' @param sensitivity (numeric, optional) Search sensitivity (`-s`
#'   parameter). Higher values are slower but more sensitive (range 1--7).
#'   If `NULL`, MMseqs2 uses its default.
#' @param min_seq_id (numeric, optional) Minimum sequence identity
#'   (0--1). If `NULL`, MMseqs2 uses its default.
#' @param e_value (numeric, optional) Maximum E-value threshold (`-e`).
#'   If `NULL`, MMseqs2 uses its default.
#' @param max_accept (integer, optional) Maximum number of hits accepted per
#'   query (`--max-accept`). Useful with `lca_mode = 1` or `4` to widen the
#'   hit set used for LCA (default: `5`).
#' @param nproc (integer) Number of threads (default: 1).
#' @param clean_pq (logical) Clean the phyloseq object before
#'   searching? (default: `TRUE`).
#' @param simplify_taxo (logical) Apply [simplify_taxo()] to the
#'   result? Only used when `behavior = "add_to_phyloseq"` (default:
#'   `TRUE`).
#' @param keep_temporary_files (logical) Keep intermediate files
#'   for debugging? (default: `FALSE`).
#' @param verbose (logical) Print progress messages? (default: `FALSE`).
#' @param cmd_args (character) Additional arguments appended to the
#'   MMseqs2 command.
#'
#' @returns
#' - If `behavior == "return_matrix"`: a [tibble][tibble::tibble] with
#'   columns `taxa_names` and one column per rank.
#' - If `behavior == "add_to_phyloseq"`: a new phyloseq object with
#'   amended `tax_table`.
#' @export
#' @author Adrien Taudière
#' @seealso [assign_blastn()], [assign_sintax()], [assign_vsearch_lca()]
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite [MMseqs2](https://github.com/soedinglab/MMseqs2):
#'   Mirdita M, Steinegger M, Breitwieser F, Soding J, Levy Karin E:
#'   Fast and sensitive taxonomic assignment to metagenomic contigs.
#'   Bioinformatics (2021).
#' @examples
#' \dontrun{
#' ref_fasta <- Biostrings::readDNAStringSet(system.file("extdata",
#'   "mini_UNITE_fungi.fasta.gz",
#'   package = "MiscMetabar", mustWork = TRUE
#' ))
#'
#' # Preferred usage: provide a SINTAX-formatted FASTA file.
#' # The function searches with easy-search and parses SINTAX headers.
#' res <- assign_mmseqs2(data_fungi_mini, ref_fasta = ref_fasta)
#' head(res)
#'
#' # Add taxonomy to phyloseq:
#' physeq_new <- assign_mmseqs2(
#'   data_fungi_mini,
#'   ref_fasta = ref_fasta,
#'   behavior = "add_to_phyloseq"
#' )
#' }
assign_mmseqs2 <- function(
  physeq = NULL,
  ref_fasta = NULL,
  database = NULL,
  seq2search = NULL,
  mmseqs2path = find_mmseqs2(),
  behavior = c("return_matrix", "add_to_phyloseq"),
  suffix = "_mmseqs2",
  lca_mode = 3,
  lca_ranks = c(
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  ),
  column_names = c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ),
  search_type = 3,
  sensitivity = NULL,
  min_seq_id = NULL,
  e_value = NULL,
  max_accept = 5,
  nproc = 1,
  clean_pq = TRUE,
  simplify_taxo = TRUE,
  keep_temporary_files = FALSE,
  verbose = FALSE,
  cmd_args = ""
) {
  behavior <- match.arg(behavior)

  if (length(lca_ranks) != length(column_names)) {
    stop("`lca_ranks` and `column_names` must have the same length.")
  }

  if (is.null(ref_fasta) && is.null(database)) {
    stop(
      "You must provide either `ref_fasta` or `database`."
    )
  }
  if (!is.null(ref_fasta) && !is.null(database)) {
    stop(
      "You set both `ref_fasta` and `database`. Please use only one."
    )
  }

  # Write query sequences to temp FASTA ----------------------------------------
  temporary_fasta_file <- file.path(tempdir(), "mmseqs2_query.fasta")
  write_temp_fasta(
    physeq = physeq,
    seq2search = seq2search,
    temporary_fasta_file = temporary_fasta_file,
    behavior = behavior,
    clean_pq = clean_pq,
    verbose = verbose
  )

  tmp_dir <- file.path(tempdir(), "mmseqs2_tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  rank_col_names <- paste0(column_names, suffix)

  # =========================================================================
  # ref_fasta path: build SINTAX taxonomy DB, then use easy-taxonomy
  # =========================================================================
  ref_fasta_path <- NULL
  if (!is.null(ref_fasta)) {
    ref_seqs <- ref_fasta
    if (is.character(ref_fasta)) {
      ref_seqs <- Biostrings::readDNAStringSet(ref_fasta)
    }

    headers <- names(ref_seqs)
    split_headers <- strsplit(headers, ";tax=", fixed = TRUE)
    ref_ids <- vapply(split_headers, \(x) trimws(x[[1]]), character(1))
    tax_strings <- vapply(
      split_headers,
      \(x) if (length(x) > 1) x[[2]] else "",
      character(1)
    )

    # Map SINTAX single-letter prefixes to rank names matching lca_ranks
    prefix_to_rank <- c(
      d = lca_ranks[1],
      k = lca_ranks[1],
      p = lca_ranks[min(2L, length(lca_ranks))],
      c = lca_ranks[min(3L, length(lca_ranks))],
      o = lca_ranks[min(4L, length(lca_ranks))],
      f = lca_ranks[min(5L, length(lca_ranks))],
      g = lca_ranks[min(6L, length(lca_ranks))],
      s = lca_ranks[min(7L, length(lca_ranks))]
    )

    # Write cleaned ref FASTA (header = ref_id only, no ;tax= suffix)
    ref_fasta_path <- file.path(tempdir(), "mmseqs2_ref.fasta")
    names(ref_seqs) <- ref_ids
    Biostrings::writeXStringSet(ref_seqs, ref_fasta_path)

    # Build MMseqs2 taxonomy database from SINTAX headers
    active_db <- .build_sintax_mmseqs2_taxdb(
      ref_ids,
      tax_strings,
      prefix_to_rank,
      ref_fasta_path,
      mmseqs2path,
      verbose = verbose
    )
  } else {
    # =========================================================================
    # database path: pre-built MMseqs2 database with NCBI taxonomy
    # =========================================================================
    active_db <- database
  }

  # Shared: run easy-taxonomy -------------------------------------------------
  output_prefix <- file.path(tempdir(), "mmseqs2_result")
  cmd <- paste(
    "easy-taxonomy",
    temporary_fasta_file,
    active_db,
    output_prefix,
    tmp_dir,
    "--lca-mode",
    lca_mode,
    "--lca-ranks",
    paste(lca_ranks, collapse = ","),
    "--search-type",
    search_type,
    "--threads",
    nproc
  )
  if (!is.null(sensitivity)) {
    cmd <- paste(cmd, "-s", sensitivity)
  }
  if (!is.null(min_seq_id)) {
    cmd <- paste(cmd, "--min-seq-id", min_seq_id)
  }
  if (!is.null(e_value)) {
    cmd <- paste(cmd, "-e", e_value)
  }
  if (!is.null(max_accept)) {
    cmd <- paste(cmd, "--max-accept", max_accept)
  }
  if (nzchar(cmd_args)) {
    cmd <- paste(cmd, cmd_args)
  }

  if (lca_mode %in% c(1, 4)) {
    message(
      "Using lca_mode = ",
      lca_mode,
      " with max_accept = ",
      max_accept,
      ". This may be slow for large datasets."
    )
  }
  if (verbose) {
    message("Running: ", mmseqs2path, " ", cmd)
  }

  mmseqs2_output <- system2(
    mmseqs2path,
    args = cmd,
    stdout = if (verbose) "" else TRUE,
    stderr = if (verbose) "" else TRUE
  )
  mmseqs2_status <- attr(mmseqs2_output, "status")
  if (!is.null(mmseqs2_status) && mmseqs2_status != 0) {
    stop(
      "MMseqs2 easy-taxonomy failed with status ",
      mmseqs2_status,
      ".\n",
      paste(mmseqs2_output, collapse = "\n")
    )
  }

  # Parse LCA output ----------------------------------------------------------
  lca_file <- paste0(output_prefix, "_lca.tsv")
  if (!file.exists(lca_file) || file.info(lca_file)$size == 0) {
    warning("No taxonomic assignments were made.")
    if (!keep_temporary_files) {
      cleanup_mmseqs2_files(temporary_fasta_file, output_prefix, tmp_dir)
    }
    if (behavior == "add_to_phyloseq") {
      return(physeq)
    } else {
      return(NULL)
    }
  }

  lca_raw <- utils::read.delim(
    lca_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    fill = TRUE,
    quote = ""
  )

  if (verbose) {
    message(
      "LCA output: ",
      ncol(lca_raw),
      " columns, ",
      nrow(lca_raw),
      " rows"
    )
  }

  lineage_col_idx <- NULL
  for (j in rev(seq_len(ncol(lca_raw)))) {
    vals <- lca_raw[[j]][!is.na(lca_raw[[j]]) & nzchar(lca_raw[[j]])]
    if (length(vals) > 0 && all(grepl(";", vals, fixed = TRUE))) {
      n_sep <- lengths(regmatches(vals[1], gregexpr(";", vals[1])))
      if (n_sep == length(lca_ranks) - 1) {
        lineage_col_idx <- j
        break
      }
    }
  }

  if (!is.null(lineage_col_idx)) {
    lineage_split <- strsplit(lca_raw[[lineage_col_idx]], ";", fixed = TRUE)
    lineage_mat <- do.call(rbind, lineage_split)
    colnames(lineage_mat) <- rank_col_names
    res_taxo <- data.frame(
      taxa_names = lca_raw[[1]],
      lineage_mat,
      stringsAsFactors = FALSE
    )
  } else {
    n_ranks <- length(lca_ranks)
    n_base <- ncol(lca_raw) - n_ranks
    if (n_base >= 1 && ncol(lca_raw) >= n_ranks + 1) {
      rank_data <- lca_raw[, (n_base + 1):ncol(lca_raw), drop = FALSE]
      colnames(rank_data) <- rank_col_names
      res_taxo <- data.frame(
        taxa_names = lca_raw[[1]],
        rank_data,
        stringsAsFactors = FALSE
      )
    } else {
      warning(
        "Could not parse MMseqs2 LCA output (",
        ncol(lca_raw),
        " columns)"
      )
      if (behavior == "add_to_phyloseq") {
        return(physeq)
      } else {
        return(NULL)
      }
    }
  }

  res_taxo <- res_taxo |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(rank_col_names),
      \(x) {
        ifelse(
          nzchar(x) &
            !x %in% c("unclassified", "unknown") &
            !startsWith(x, "uc_"),
          x,
          NA_character_
        )
      }
    ))

  # Cleanup -------------------------------------------------------------------
  if (keep_temporary_files) {
    message("Temporary files are located at ", tempdir())
  } else {
    if (!is.null(ref_fasta_path)) {
      unlink(ref_fasta_path)
      unlink(
        list.files(tempdir(), pattern = "^mmseqs2_refDB", full.names = TRUE),
        recursive = TRUE
      )
      unlink(file.path(tempdir(), "mmseqs2_taxdump"), recursive = TRUE)
      unlink(file.path(tempdir(), "mmseqs2_seqtaxid.tsv"))
      unlink(
        file.path(tempdir(), "mmseqs2_createtaxdb_tmp"),
        recursive = TRUE
      )
    }
    cleanup_mmseqs2_files(temporary_fasta_file, output_prefix, tmp_dir)
  }

  # Return results (shared) ---------------------------------------------------
  if (behavior == "return_matrix") {
    return(res_taxo)
  } else if (behavior == "add_to_phyloseq") {
    tax_tab <- as.data.frame(as.matrix(physeq@tax_table))
    tax_tab$taxa_names <- taxa_names(physeq)

    new_physeq <- physeq

    new_tax_tab <- dplyr::left_join(
      tax_tab,
      res_taxo,
      by = dplyr::join_by(taxa_names)
    ) |>
      dplyr::select(-taxa_names) |>
      as.matrix()

    new_physeq@tax_table <- tax_table(new_tax_tab)
    taxa_names(new_physeq@tax_table) <- taxa_names(physeq)

    if (simplify_taxo) {
      new_physeq <- simplify_taxo(new_physeq)
    }

    verify_pq(new_physeq)

    return(new_physeq)
  }
}

################################################################################
#' Recluster sequences of a phyloseq object or cluster a list of DNA sequences
#'   using MMseqs2 software
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' A wrapper of the [MMseqs2](https://github.com/soedinglab/MMseqs2)
#' `easy-cluster` or `easy-linclust` workflow.
#'
#' @inheritParams clean_pq
#' @param dna_seq You may directly use a character vector of DNA sequences
#'   in place of physeq args. When physeq is set, dna sequences take the
#'   value of `physeq@refseq`.
#' @param nproc (default: 1) Number of threads.
#' @param id (default: 0.97) Minimum sequence identity threshold (0--1).
#' @param mmseqs2path Path to the `mmseqs` binary
#'   (default: [find_mmseqs2()]).
#' @param tax_adjust (Default 0) See the man page
#'   of [merge_taxa_vec()] for more details.
#'   To conserve the taxonomic rank of the most abundant taxa (ASV, OTU, ...),
#'   set `tax_adjust` to 0 (default).
#' @param rank_propagation (logical, default FALSE). Do we propagate the
#'   NA value from lower taxonomic rank to upper rank?
#'   See the man page of [merge_taxa_vec()] for more details.
#' @param mmseqs2_cluster_method (default: `"easy-cluster"`) Either
#'   `"easy-cluster"` (cascaded clustering, more sensitive) or
#'   `"easy-linclust"` (linear-time clustering, faster for huge datasets).
#' @param coverage (numeric, default: 0.8) Alignment coverage threshold
#'   (0--1), passed to `-c`.
#' @param cov_mode (integer, default: 0) Coverage mode:
#'   - `0`: `alnRes / max(qLen, tLen)`
#'   - `1`: `alnRes / tLen`
#'   - `2`: `alnRes / qLen`
#' @param cluster_mode (integer, default: 0) Clustering algorithm:
#'   - `0`: greedy set cover (default)
#'   - `1`: connected components
#'   - `2`: greedy incremental
#' @param mmseqs2_args (character, default: `""`) Additional arguments
#'   passed to the MMseqs2 clustering command.
#' @param keep_temporary_files (logical, default: FALSE) Keep intermediate
#'   files for debugging?
#'
#' @return A new object of class `physeq` or a data.frame of cluster
#'   membership if `dna_seq` was used.
#' @seealso [postcluster_pq()], [vsearch_clustering()], [swarm_clustering()]
#' @export
#' @author Adrien Taudière
#' @references
#'   MMseqs2 can be downloaded from
#'   \url{https://github.com/soedinglab/MMseqs2}.
#'   More information in the associated publication
#'   \doi{10.1038/nbt.3988}.
#' @examplesIf MiscMetabar::is_mmseqs2_installed()
#' \donttest{
#' d_mm <- mmseqs2_clustering(data_fungi_mini)
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite [MMseqs2](https://github.com/soedinglab/MMseqs2).
mmseqs2_clustering <- function(
  physeq = NULL,
  dna_seq = NULL,
  nproc = 1,
  id = 0.97,
  mmseqs2path = find_mmseqs2(),
  tax_adjust = 0,
  rank_propagation = FALSE,
  mmseqs2_cluster_method = c("easy-cluster", "easy-linclust"),
  coverage = 0.8,
  cov_mode = 0,
  cluster_mode = 0,
  mmseqs2_args = "",
  keep_temporary_files = FALSE
) {
  mmseqs2_cluster_method <- match.arg(mmseqs2_cluster_method)
  dna <- physeq_or_string_to_dna(physeq = physeq, dna_seq = dna_seq)

  temporary_fasta_file <- paste0(tempdir(), "/mmseqs2_clust_input.fasta")
  Biostrings::writeXStringSet(dna, temporary_fasta_file)

  output_prefix <- paste0(tempdir(), "/mmseqs2_clust")
  tmp_dir <- paste0(tempdir(), "/mmseqs2_clust_tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  # Build command
  cmd <- paste0(
    " ",
    mmseqs2_cluster_method,
    " ",
    temporary_fasta_file,
    " ",
    output_prefix,
    " ",
    tmp_dir,
    " --min-seq-id ",
    id,
    " -c ",
    coverage,
    " --cov-mode ",
    cov_mode,
    " --cluster-mode ",
    cluster_mode,
    " --threads ",
    nproc
  )

  if (nzchar(mmseqs2_args)) {
    cmd <- paste0(cmd, " ", mmseqs2_args)
  }

  mmseqs2_output <- system2(
    mmseqs2path,
    args = cmd,
    stdout = TRUE,
    stderr = TRUE
  )
  mmseqs2_status <- attr(mmseqs2_output, "status")

  if (!is.null(mmseqs2_status) && mmseqs2_status != 0) {
    stop(
      "MMseqs2 ",
      mmseqs2_cluster_method,
      " failed with status ",
      mmseqs2_status,
      ".\n",
      paste(mmseqs2_output, collapse = "\n")
    )
  }

  # Parse cluster TSV: columns are representative_id and member_id
  cluster_tsv <- paste0(output_prefix, "_cluster.tsv")
  if (!file.exists(cluster_tsv) || file.info(cluster_tsv)$size == 0) {
    stop("MMseqs2 clustering produced no output.")
  }

  clust_tab <- utils::read.delim(
    cluster_tsv,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("representative", "member")
  )

  # Convert to numeric cluster vector (same format as vsearch_clustering)
  rep_ids <- unique(clust_tab$representative)
  cluster_map <- stats::setNames(seq_along(rep_ids), rep_ids)
  clusters <- cluster_map[clust_tab$representative]
  names(clusters) <- clust_tab$member

  # Cleanup
  if (keep_temporary_files) {
    message("Temporary files are located at ", tempdir())
  } else {
    unlink(temporary_fasta_file)
    cleanup_mmseqs2_files(temporary_fasta_file, output_prefix, tmp_dir)
  }

  if (inherits(physeq, "phyloseq")) {
    clusters <- clusters[match(taxa_names(physeq), names(clusters))]
    new_obj <- merge_taxa_vec(
      physeq,
      clusters,
      tax_adjust = tax_adjust,
      rank_propagation = rank_propagation
    )
  } else if (inherits(dna_seq, "character")) {
    new_obj <- clust_tab
  } else {
    stop(
      "You must set the args physeq (object of class phyloseq) ",
      "or dna_seq (character vector)."
    )
  }

  return(new_obj)
}

################################################################################
#' Remove temporary MMseqs2 files
#' @param fasta_file Path to the temporary FASTA query file.
#' @param output_prefix Prefix used for MMseqs2 output files.
#' @param tmp_dir Path to the MMseqs2 temporary directory.
#' @returns NULL (called for side effects).
#' @noRd
cleanup_mmseqs2_files <- function(fasta_file, output_prefix, tmp_dir) {
  unlink(fasta_file)
  output_files <- list.files(
    dirname(output_prefix),
    pattern = paste0("^", basename(output_prefix)),
    full.names = TRUE
  )
  unlink(output_files)
  unlink(tmp_dir, recursive = TRUE)
}
################################################################################
