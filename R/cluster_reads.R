################################################################################
#' Cluster raw reads into OTUs (vsearch or SWARM), no DADA2 required
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Cluster marker-gene reads into similarity-threshold OTUs and assemble
#' them into a [phyloseq::phyloseq-class()] object with an `otu_table` of
#' per-sample abundances and the representative (seed/centroid) sequence of
#' each OTU in the `refseq` slot: reads are dereplicated per sample with
#' [vsearch](https://github.com/torognes/vsearch), pooled, dereplicated
#' again, and clustered with vsearch (`--cluster_size`, default) or
#' [SWARM](https://github.com/torognes/swarm) (`d > 1` recommended).
#'
#' This complements [asv2otu()], which post-clusters ASVs already stored in
#' a phyloseq object: `cluster_reads()` starts one step earlier, directly
#' from sequence files or a vector of DNA sequences. Use [denoised_reads()]
#' instead of this function to denoise reads into ASVs/zOTUs (dada2
#' pipeline, SWARM `d = 1` or vsearch UNOISE).
#'
#' @param path_to_fastx (character) Path(s) to fasta or fastq files
#'   (optionally gzipped). Each file is treated as one sample; the sample
#'   name is the basename stripped of its extension. Set either
#'   `path_to_fastx` or `dna_seq`, not both.
#' @param dna_seq A named character vector or
#'   [Biostrings::DNAStringSet-class()] of DNA sequences, treated as one
#'   sample. Duplicated sequences are counted; alternatively, already
#'   dereplicated sequences may carry their abundance as a `;size=N` suffix
#'   in their names. Set either `path_to_fastx` or `dna_seq`, not both.
#' @param sample_name (default: `"sample_1"`) Sample name used when
#'   `dna_seq` is set.
#' @param method (default: `"vsearch"`) Clustering tool: `"vsearch"` or
#'   `"swarm"`.
#' @param id (default: 0.97) Identity threshold for vsearch clustering.
#' @param d (default: 1) SWARM resolution. `d = 1` is a denoising
#'   resolution and triggers a warning pointing to [denoised_reads()]; use
#'   `d > 1` for OTU-like clustering.
#' @param fastidious (logical, default: TRUE) Use SWARM fastidious mode.
#'   Ignored when `d != 1`.
#' @param minsize (default: 1) Discard dereplicated sequences with an
#'   abundance strictly lower than `minsize` before clustering. The common
#'   pre-clustering filter `minsize = 2` removes singletons.
#' @param nproc (default: 1) Number of threads for the clustering tool.
#' @param vsearchpath (default: [find_vsearch()]) Path to the vsearch binary.
#' @param swarmpath (default: `"swarm"`) Path to the swarm binary.
#' @param vsearch_args (default: `"--strand both"`) Extra arguments passed to
#'   vsearch clustering.
#' @param swarm_args (default: `""`) Extra arguments passed to swarm.
#' @param keep_temporary_files (logical, default: FALSE) Keep the
#'   dereplicated fasta files, the pooled fasta, the cluster table and the
#'   representative-sequences fasta in a temporary directory.
#'
#' @returns A [phyloseq::phyloseq-class()] object with one row (taxon) per
#'   OTU named `OTU_1`, `OTU_2`, ... ordered by decreasing total abundance,
#'   one column (sample) per input file, the representative sequence of each
#'   OTU in the `refseq` slot, and (for file input) the source file paths in
#'   the `sample_data` slot.
#'
#' @details
#' vsearch is required (per-sample and pool-level dereplication, and
#' clustering for `method = "vsearch"`); swarm is additionally required for
#' `method = "swarm"`. Use [is_vsearch_installed()] and
#' [is_swarm_installed()] to check availability. This function is mainly a
#' wrapper of the work of others: please cite
#' [vsearch](https://github.com/torognes/vsearch)
#' (\url{https://pubmed.ncbi.nlm.nih.gov/27781170}) and/or
#' [SWARM](https://github.com/torognes/swarm)
#' (\doi{10.1093/bioinformatics/btab493}).
#'
#' @export
#' @author Adrien Taudière
#' @seealso [denoised_reads()], [asv2otu()], [swarm_clustering()],
#'   [vsearch_clustering()], [find_vsearch()], [install_vsearch()]
#'
#' @examplesIf MiscMetabar::is_vsearch_installed() && MiscMetabar::is_swarm_installed()
#' f1 <- system.file("extdata", "ex_R1_001.fasta", package = "MiscMetabar")
#' f2 <- system.file("extdata", "ex_R1_002.fasta", package = "MiscMetabar")
#' (pq_otu <- cluster_reads(path_to_fastx = c(f1, f2)))
#'
#' cluster_reads(path_to_fastx = c(f1, f2), method = "swarm", d = 3)
cluster_reads <- function(
  path_to_fastx = NULL,
  dna_seq = NULL,
  sample_name = "sample_1",
  method = c("vsearch", "swarm"),
  id = 0.97,
  d = 1,
  fastidious = TRUE,
  minsize = 1,
  nproc = 1,
  vsearchpath = find_vsearch(),
  swarmpath = "swarm",
  vsearch_args = "--strand both",
  swarm_args = "",
  keep_temporary_files = FALSE
) {
  method <- match.arg(method)
  sample_names <- reads_validate_inputs(path_to_fastx, dna_seq, sample_name)

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
  if (method == "swarm" && d == 1) {
    cli::cli_warn(
      "SWARM with {.code d = 1} is a denoising resolution (zOTU-like \\
      output); consider {.fn MiscMetabar::denoised_reads} instead."
    )
  }

  tmp <- tempfile("cluster_reads_")
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

  if (method == "swarm") {
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
  } else {
    run_tool(
      vsearchpath,
      paste0(
        " --cluster_size ",
        shQuote(pool$unique_file),
        " ",
        vsearch_args,
        " -id ",
        id,
        " --sizein --minseqlength 1 --centroids ",
        shQuote(rep_file),
        " --uc ",
        shQuote(uc_file)
      ),
      "vsearch"
    )
  }

  if (!file.exists(uc_file)) {
    cli::cli_abort(
      "Clustering failed: no cluster table was produced. \\
      Set {.code keep_temporary_files = TRUE} to inspect intermediate files."
    )
  }

  new_physeq <- reads_build_phyloseq(
    uc_file,
    pool$map_file,
    rep_file,
    sample_names,
    path_to_fastx,
    taxa_prefix = "OTU_"
  )

  if (!keep_temporary_files) {
    unlink(tmp, recursive = TRUE)
  } else {
    cli::cli_inform("Temporary files kept in {.file {tmp}}")
  }

  new_physeq
}
################################################################################
