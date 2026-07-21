################################################################################
# Internal pipeline helpers shared by denoised_reads() and cluster_reads().
# Not exported.
################################################################################

# Run an external tool and abort with its last output lines on failure.
run_tool <- function(cmd, args, tool_name) {
  out <- system2(cmd, args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cli::cli_abort(
      c(
        "{.pkg {tool_name}} failed with status {status}:",
        x = paste(utils::tail(out, 5), collapse = "\n")
      )
    )
  }
  invisible(out)
}

# Validate the shared inputs of denoised_reads() and cluster_reads() and
# derive the sample names (from file basenames, or sample_name for dna_seq).
reads_validate_inputs <- function(path_to_fastx, dna_seq, sample_name) {
  if (!xor(is.null(path_to_fastx), is.null(dna_seq))) {
    cli::cli_abort(
      "You must set exactly one of {.arg path_to_fastx} or {.arg dna_seq}."
    )
  }
  if (!is.null(path_to_fastx)) {
    missing_files <- path_to_fastx[!file.exists(path_to_fastx)]
    if (length(missing_files) > 0) {
      cli::cli_abort("File{?s} not found: {.file {missing_files}}")
    }
    sample_names <- gsub(
      "[^A-Za-z0-9_.-]",
      "_",
      sub(
        "\\.(fastq|fq|fasta|fa|fna)(\\.gz)?$",
        "",
        basename(path_to_fastx)
      )
    )
    if (anyDuplicated(sample_names)) {
      cli::cli_abort(
        "Sample names derived from file basenames are not unique: \\
        {.val {unique(sample_names[duplicated(sample_names)])}}."
      )
    }
  } else {
    sample_names <- gsub("[^A-Za-z0-9_.-]", "_", sample_name)
  }
  sample_names
}

# Write per-sample dereplicated fasta files (headers "<sample>_<i>;size=N")
# into tmp and return their paths: vsearch dereplication for file input,
# R-side counting for dna_seq (with optional ";size=N" name annotations).
reads_dereplicate_samples <- function(
  vsearchpath,
  path_to_fastx,
  dna_seq,
  sample_names,
  minsize,
  tmp
) {
  if (!is.null(path_to_fastx)) {
    derep_files <- file.path(tmp, paste0("derep_", sample_names, ".fasta"))
    for (i in seq_along(path_to_fastx)) {
      run_tool(
        vsearchpath,
        paste0(
          " --derep_fulllength ",
          shQuote(path_to_fastx[i]),
          " --sizein --sizeout --fasta_width 0 --minseqlength 1",
          " --minuniquesize ",
          minsize,
          " --relabel ",
          sample_names[i],
          "_ --output ",
          shQuote(derep_files[i])
        ),
        "vsearch"
      )
    }
    return(derep_files)
  }

  seqs <- as.character(dna_seq)
  seq_names <- names(dna_seq)
  has_size <- !is.null(seq_names) &&
    all(grepl(";size=[0-9]+$", seq_names)) &&
    !anyDuplicated(seqs)
  if (has_size) {
    abund <- as.integer(sub("^.*;size=", "", seq_names))
  } else {
    tbl <- table(seqs)
    seqs <- names(tbl)
    abund <- as.integer(tbl)
  }
  keep <- abund >= minsize
  seqs <- seqs[keep]
  abund <- abund[keep]
  if (length(seqs) == 0) {
    cli::cli_abort(
      "No sequence left after the {.code minsize = {minsize}} filter."
    )
  }
  ord <- order(abund, decreasing = TRUE)
  dna <- Biostrings::DNAStringSet(seqs[ord])
  names(dna) <- paste0(
    sample_names,
    "_",
    seq_along(dna),
    ";size=",
    abund[ord]
  )
  derep_file <- file.path(tmp, paste0("derep_", sample_names, ".fasta"))
  Biostrings::writeXStringSet(dna, derep_file, width = 20000)
  derep_file
}

# The same sequence can occur in several samples, so the pool of per-sample
# dereplicated files is not itself dereplicated: swarm and vsearch
# clustering require a fully dereplicated input. A second, pool-level
# dereplication (with --uc) both produces that input and records, for every
# pooled entry, which representative sequence it was merged into -- this is
# what allows per-sample counts to be rebuilt.
reads_pool_dereplicate <- function(vsearchpath, derep_files, tmp) {
  pool_file <- file.path(tmp, "pool.fasta")
  file.create(pool_file)
  for (f in derep_files) {
    file.append(pool_file, f)
  }
  unique_file <- file.path(tmp, "unique.fasta")
  map_file <- file.path(tmp, "derep_map.uc")
  run_tool(
    vsearchpath,
    paste0(
      " --derep_fulllength ",
      shQuote(pool_file),
      " --sizein --sizeout --fasta_width 0 --minseqlength 1",
      " --uc ",
      shQuote(map_file),
      " --output ",
      shQuote(unique_file)
    ),
    "vsearch"
  )
  list(unique_file = unique_file, map_file = map_file)
}

# Parse a cluster table (uclust format), the dereplication map and the
# representative-sequences fasta, and assemble the phyloseq object.
reads_build_phyloseq <- function(
  uc_file,
  map_file,
  rep_file,
  sample_names,
  path_to_fastx,
  taxa_prefix
) {
  read_uc <- function(path) {
    uc_lines <- strsplit(readLines(path), "\t", fixed = TRUE)
    uc <- do.call(
      rbind,
      lapply(uc_lines, function(x) {
        length(x) <- 10
        x
      })
    )
    uc <- as.data.frame(uc, stringsAsFactors = FALSE)
    colnames(uc) <- c(
      "type",
      "cluster",
      "width",
      "identity",
      "strand",
      "f6",
      "f7",
      "cigar",
      "query",
      "target"
    )
    uc[uc$type %in% c("S", "H"), ]
  }
  strip_size <- function(x) {
    sub(";size=[0-9]+;?$", "", x)
  }
  get_size <- function(x) {
    size <- suppressWarnings(
      as.integer(sub("^.*;size=([0-9]+);?$", "\\1", x))
    )
    size[is.na(size)] <- 1L
    size
  }

  uc <- read_uc(uc_file)
  if (nrow(uc) == 0) {
    cli::cli_abort("Clustering produced an empty cluster table.")
  }
  cluster_of_unique <- as.integer(uc$cluster)
  names(cluster_of_unique) <- strip_size(uc$query)

  map <- read_uc(map_file)
  map_unique <- ifelse(map$type == "S", map$query, map$target)
  map_unique <- strip_size(map_unique)
  map_sample <- strip_size(map$query)
  map_size <- get_size(map$query)

  sample_of <- rep(NA_character_, length(map_sample))
  for (s in sample_names[order(nchar(sample_names), decreasing = TRUE)]) {
    hit <- startsWith(map_sample, paste0(s, "_"))
    sample_of[hit & is.na(sample_of)] <- s
  }
  if (any(is.na(sample_of))) {
    cli::cli_abort(
      "Could not map cluster records back to samples. \\
      Set {.code keep_temporary_files = TRUE} to inspect intermediate files."
    )
  }

  map_cluster <- cluster_of_unique[map_unique]
  if (any(is.na(map_cluster))) {
    cli::cli_abort(
      "Could not map every sequence to a cluster. \\
      Set {.code keep_temporary_files = TRUE} to inspect intermediate files."
    )
  }

  mat <- tapply(
    map_size,
    list(
      cluster = factor(map_cluster),
      sample = factor(sample_of, levels = sample_names)
    ),
    sum,
    default = 0
  )
  mat <- unclass(mat)
  mat[is.na(mat)] <- 0
  storage.mode(mat) <- "integer"

  otu_order <- order(rowSums(mat), decreasing = TRUE)
  cluster_ids <- as.integer(rownames(mat))[otu_order]
  taxa_names <- paste0(taxa_prefix, seq_along(otu_order))
  mat <- mat[otu_order, , drop = FALSE]
  rownames(mat) <- taxa_names

  is_seed <- uc$type == "S"
  seed_of_cluster <- strip_size(uc$query[is_seed])[
    match(cluster_ids, as.integer(uc$cluster[is_seed]))
  ]
  rep_seqs <- Biostrings::readDNAStringSet(rep_file)
  rep_key <- strip_size(names(rep_seqs))
  seq_idx <- match(seed_of_cluster, rep_key)
  if (any(is.na(seq_idx))) {
    cli::cli_abort(
      "Could not retrieve a representative sequence for every taxon. \\
      Set {.code keep_temporary_files = TRUE} to inspect intermediate files."
    )
  }
  refseq <- rep_seqs[seq_idx]
  names(refseq) <- taxa_names

  if (!is.null(path_to_fastx)) {
    sam <- data.frame(
      source_file = path_to_fastx,
      row.names = sample_names
    )
    phyloseq::phyloseq(
      phyloseq::otu_table(mat, taxa_are_rows = TRUE),
      phyloseq::sample_data(sam),
      refseq
    )
  } else {
    phyloseq::phyloseq(
      phyloseq::otu_table(mat, taxa_are_rows = TRUE),
      refseq
    )
  }
}
################################################################################

################################################################################
#' Denoise raw reads into ASVs / zOTUs (dada2, SWARM or vsearch UNOISE)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Denoise marker-gene reads into amplicon sequence variants (ASVs, named
#' `ASV_1`, `ASV_2`, ...) and assemble them into a
#' [phyloseq::phyloseq-class()] object with an `otu_table` of per-sample
#' abundances and the representative sequence of each variant in the
#' `refseq` slot. Three denoising engines are available:
#'
#' - `method = "vsearch"` (default): [vsearch](https://github.com/torognes/vsearch)
#'   UNOISE3 denoising (`--cluster_unoise`), producing zOTUs. Works on
#'   fasta, fastq and `dna_seq` inputs.
#' - `method = "swarm"`: [SWARM](https://github.com/torognes/swarm) with
#'   `d = 1` (fastidious by default). Works on fasta, fastq and `dna_seq`
#'   inputs.
#' - `method = "dada2"`: the classical [dada2](https://benjjneb.github.io/dada2/)
#'   denoising pipeline (`dada2::dada()` with `selfConsist = TRUE` by
#'   default, overridable through `...`), followed by de-novo chimera
#'   removal (`dada2::removeBimeraDenovo()`). Requires fastq files (quality
#'   scores drive the dada2 error model); reads containing non-ACGT bases
#'   are discarded with a message.
#'
#' For fasta/fastq input, reads are first dereplicated per sample with
#' vsearch, then pooled and dereplicated again so that abundances per
#' sample can be rebuilt after denoising (see the details of
#' [cluster_reads()]). Use [cluster_reads()] instead of this function to
#' cluster reads into similarity-threshold OTUs (vsearch `--cluster_size`
#' or SWARM with `d > 1`).
#'
#' @param path_to_fastx (character) Path(s) to fasta or fastq files
#'   (optionally gzipped). Each file is treated as one sample; the sample
#'   name is the basename stripped of its extension. Set either
#'   `path_to_fastx` or `dna_seq`, not both.
#' @param dna_seq A named character vector or
#'   [Biostrings::DNAStringSet-class()] of DNA sequences, treated as one
#'   sample. Duplicated sequences are counted; alternatively, already
#'   dereplicated sequences may carry their abundance as a `;size=N` suffix
#'   in their names. Not available for `method = "dada2"` (no quality
#'   scores). Set either `path_to_fastx` or `dna_seq`, not both.
#' @param sample_name (default: `"sample_1"`) Sample name used when
#'   `dna_seq` is set.
#' @param method (default: `"vsearch"`) Denoising engine: `"vsearch"`
#'   (UNOISE3), `"swarm"` or `"dada2"`.
#' @param unoise_alpha (numeric, default 2) The `alpha` parameter of
#'   vsearch UNOISE (trade-off between sensitivity and specificity).
#' @param d (default: 1) SWARM resolution. `d = 1` is the denoising
#'   resolution; larger values cluster at a lower resolution and trigger a
#'   warning pointing to [cluster_reads()].
#' @param fastidious (logical, default: TRUE) Use SWARM fastidious mode.
#'   Ignored when `d != 1`.
#' @param minsize (default: 1) Discard dereplicated sequences with an
#'   abundance strictly lower than `minsize` before denoising. Also passed
#'   to vsearch UNOISE (`--minsize`).
#' @param nproc (default: 1) Number of threads (SWARM and dada2 chimera
#'   removal).
#' @param remove_chimeras (logical, default: TRUE) For `method = "dada2"`,
#'   remove chimeras de novo with `dada2::removeBimeraDenovo()` after
#'   denoising.
#' @param vsearchpath (default: [find_vsearch()]) Path to the vsearch binary.
#' @param swarmpath (default: `"swarm"`) Path to the swarm binary.
#' @param swarm_args (default: `""`) Extra arguments passed to swarm.
#' @param keep_temporary_files (logical, default: FALSE) Keep the
#'   dereplicated fasta files, the pooled fasta, the cluster table and the
#'   representative-sequences fasta in a temporary directory.
#' @param ... Additional arguments passed to [dada2::dada()] when
#'   `method = "dada2"` (e.g. `err`, `selfConsist`, `OMEGA_A`, `pool`).
#'
#' @returns A [phyloseq::phyloseq-class()] object with one row (taxon) per
#'   ASV/zOTU named `ASV_1`, `ASV_2`, ... ordered by decreasing total
#'   abundance, one column (sample) per input file, the representative
#'   sequence of each variant in the `refseq` slot, and (for file input)
#'   the source file paths in the `sample_data` slot.
#'
#' @details
#' vsearch is required for all methods (per-sample dereplication and, for
#' `method = "vsearch"`, UNOISE denoising); swarm is additionally required
#' for `method = "swarm"`. Use [is_vsearch_installed()] and
#' [is_swarm_installed()] to check availability. This function is mainly a
#' wrapper of the work of others: please cite
#' [vsearch](https://github.com/torognes/vsearch)
#' (\url{https://pubmed.ncbi.nlm.nih.gov/27781170}),
#' [SWARM](https://github.com/torognes/swarm)
#' (\doi{10.1093/bioinformatics/btab493}) and/or
#' [dada2](https://doi.org/10.1038/nmeth.3869) depending on the method used.
#'
#' @export
#' @author Adrien Taudière
#' @seealso [cluster_reads()], [asv2otu()], [swarm_clustering()],
#'   [vsearch_clustering()], [chimera_removal_vs()]
#'
#' @examplesIf MiscMetabar::is_vsearch_installed() && MiscMetabar::is_swarm_installed()
#' f1 <- system.file("extdata", "ex_R1_001.fasta", package = "MiscMetabar")
#' f2 <- system.file("extdata", "ex_R1_002.fasta", package = "MiscMetabar")
#' (pq_zotu <- denoised_reads(path_to_fastx = c(f1, f2)))
#'
#' denoised_reads(path_to_fastx = c(f1, f2), method = "swarm")
#' @examples
#' \donttest{
#' if (requireNamespace("dada2")) {
#'   fq1 <- system.file("extdata", "ex_R1_001.fastq.gz", package = "MiscMetabar")
#'   fq2 <- system.file("extdata", "ex.fastq", package = "MiscMetabar")
#'   denoised_reads(path_to_fastx = c(fq1, fq2), method = "dada2")
#' }
#' }
denoised_reads <- function(
  path_to_fastx = NULL,
  dna_seq = NULL,
  sample_name = "sample_1",
  method = c("vsearch", "swarm", "dada2"),
  unoise_alpha = 2,
  d = 1,
  fastidious = TRUE,
  minsize = 1,
  nproc = 1,
  remove_chimeras = TRUE,
  vsearchpath = find_vsearch(),
  swarmpath = "swarm",
  swarm_args = "",
  keep_temporary_files = FALSE,
  ...
) {
  method <- match.arg(method)
  sample_names <- reads_validate_inputs(path_to_fastx, dna_seq, sample_name)

  if (method == "dada2") {
    if (!is.null(dna_seq)) {
      cli::cli_abort(
        c(
          "{.code method = \"dada2\"} requires fastq files: the dada2 \\
          error model needs quality scores, which a vector of sequences \\
          does not carry.",
          i = "Use {.code method = \"vsearch\"} or \\
          {.code method = \"swarm\"} for {.arg dna_seq} input."
        )
      )
    }
    is_fastq <- grepl(
      "\\.(fastq|fq)(\\.gz)?$",
      path_to_fastx,
      ignore.case = TRUE
    )
    if (!all(is_fastq)) {
      cli::cli_abort(
        c(
          "{.code method = \"dada2\"} requires fastq files (quality \\
          scores drive the dada2 error model); not fastq: \\
          {.file {path_to_fastx[!is_fastq]}}",
          i = "Use {.code method = \"vsearch\"} or \\
          {.code method = \"swarm\"} for fasta input."
        )
      )
    }
    return(denoised_reads_dada2(
      path_to_fastx = path_to_fastx,
      sample_names = sample_names,
      nproc = nproc,
      remove_chimeras = remove_chimeras,
      ...
    ))
  }

  if (!is_vsearch_installed(vsearchpath)) {
    cli::cli_abort(
      c(
        "vsearch is required but was not found.",
        i = "Install it with {.fn MiscMetabar::install_vsearch} or set \\
        {.arg vsearchpath}."
      )
    )
  }
  if (method == "swarm" && !is_swarm_installed(swarmpath)) {
    cli::cli_abort(
      c(
        "swarm is required for {.code method = \"swarm\"} but was not found.",
        i = "Install swarm (<https://github.com/torognes/swarm>) or use \\
        {.code method = \"vsearch\"}."
      )
    )
  }
  if (method == "swarm" && d > 1) {
    cli::cli_warn(
      "SWARM with {.code d > 1} clusters below denoising resolution \\
      (OTU-like output); consider {.fn MiscMetabar::cluster_reads} instead."
    )
  }

  tmp <- tempfile("denoised_reads_")
  dir.create(tmp, showWarnings = FALSE)

  derep_files <- reads_dereplicate_samples(
    vsearchpath,
    path_to_fastx,
    dna_seq,
    sample_names,
    minsize,
    tmp
  )
  pool <- reads_pool_dereplicate(vsearchpath, derep_files, tmp)

  uc_file <- file.path(tmp, "clusters.uc")
  rep_file <- file.path(tmp, "representatives.fasta")

  if (method == "vsearch") {
    run_tool(
      vsearchpath,
      paste0(
        " --cluster_unoise ",
        shQuote(pool$unique_file),
        " --unoise_alpha ",
        unoise_alpha,
        " --minsize ",
        minsize,
        " --sizein --minseqlength 1 --centroids ",
        shQuote(rep_file),
        " --uc ",
        shQuote(uc_file)
      ),
      "vsearch"
    )
  } else {
    if (d != 1) {
      fastidious <- FALSE
    }
    if (fastidious && !grepl("--fastidious", swarm_args)) {
      swarm_args <- paste(swarm_args, "--fastidious")
    }
    run_tool(
      swarmpath,
      paste0(
        shQuote(pool$unique_file),
        " -o ",
        shQuote(file.path(tmp, "clusters.txt")),
        " -u ",
        shQuote(uc_file),
        " -w ",
        shQuote(rep_file),
        " -z -t ",
        nproc,
        " -d ",
        d,
        " ",
        swarm_args
      ),
      "swarm"
    )
  }

  if (!file.exists(uc_file)) {
    cli::cli_abort(
      "Denoising failed: no cluster table was produced. \\
      Set {.code keep_temporary_files = TRUE} to inspect intermediate files."
    )
  }

  new_physeq <- reads_build_phyloseq(
    uc_file,
    pool$map_file,
    rep_file,
    sample_names,
    path_to_fastx,
    taxa_prefix = "ASV_"
  )

  if (!keep_temporary_files) {
    unlink(tmp, recursive = TRUE)
  } else {
    cli::cli_inform("Temporary files kept in {.file {tmp}}")
  }

  new_physeq
}
################################################################################

################################################################################
# Classical dada2 denoising pipeline for denoised_reads(method = "dada2").
# Per-sample derepFastq + dada (selfConsist error estimation by default,
# overridable through ...), makeSequenceTable and optional de-novo chimera
# removal. Not exported.
denoised_reads_dada2 <- function(
  path_to_fastx,
  sample_names,
  nproc,
  remove_chimeras,
  ...
) {
  dots <- list(...)

  dd_list <- vector("list", length(path_to_fastx))
  names(dd_list) <- sample_names
  for (i in seq_along(path_to_fastx)) {
    derep <- dada2::derepFastq(path_to_fastx[i])
    not_acgt <- grepl("[^ACGT]", names(derep$uniques))
    if (any(not_acgt)) {
      cli::cli_inform(
        "{.val {sample_names[i]}}: discarding {sum(not_acgt)} unique \\
        sequence{?s} containing non-ACGT bases before dada2 denoising."
      )
      derep$uniques <- derep$uniques[!not_acgt]
      derep$quals <- derep$quals[!not_acgt, , drop = FALSE]
    }
    dada_args <- utils::modifyList(
      list(derep = derep, err = NULL, selfConsist = TRUE, verbose = FALSE),
      dots
    )
    dd_list[[i]] <- tryCatch(
      do.call(dada2::dada, dada_args),
      error = function(e) {
        cli::cli_abort(
          "{.fn dada2::dada} failed on sample \\
          {.val {sample_names[i]}}: {conditionMessage(e)}"
        )
      }
    )
  }

  seqtab <- dada2::makeSequenceTable(dd_list)
  if (remove_chimeras) {
    seqtab <- dada2::removeBimeraDenovo(
      seqtab,
      method = "consensus",
      multithread = nproc > 1,
      verbose = FALSE
    )
  }

  mat <- t(seqtab)
  keep <- rowSums(mat) > 0
  if (!any(keep)) {
    cli::cli_abort("No sequence variant left after dada2 denoising.")
  }
  mat <- mat[keep, , drop = FALSE]
  storage.mode(mat) <- "integer"

  asv_order <- order(rowSums(mat), decreasing = TRUE)
  asv_names <- paste0("ASV_", seq_len(nrow(mat)))
  refseq <- Biostrings::DNAStringSet(rownames(mat)[asv_order])
  names(refseq) <- asv_names
  mat <- mat[asv_order, , drop = FALSE]
  rownames(mat) <- asv_names

  sam <- data.frame(
    source_file = path_to_fastx,
    row.names = sample_names
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(mat, taxa_are_rows = TRUE),
    phyloseq::sample_data(sam),
    refseq
  )
}
################################################################################
