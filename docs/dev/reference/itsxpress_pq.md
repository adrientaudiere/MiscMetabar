# Trim FASTQ files to the ITS region with ITSxpress before DADA2 denoising

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Run [ITSxpress](https://github.com/USDA-ARS-GBRU/itsxpress) (Rivers et
al. 2018; Einarsson & Rivers 2024) on paired-end or single-end FASTQ
files to trim off the conserved ribosomal flanks (SSU, 5.8S, LSU) while
preserving per-base Phred quality scores. Unlike
[`itsx_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsx_pq.md),
which runs ITSx on already-denoised FASTA reference sequences in a
`phyloseq` object, `itsxpress_pq()` operates directly on raw FASTQ reads
**before** DADA2 denoising
([`filter_trim()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/filter_trim.md),
[`denoised_reads()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/denoised_reads.md)).

By default (`paired_end = TRUE` and `single_file_output = FALSE`),
`itsxpress_pq()` writes unmerged forward and reverse trimmed FASTQ files
(`--outfile` and `--outfile2`), which is the layout required for DADA2's
paired-end error model.

## Usage

``` r
itsxpress_pq(
  path_to_fastq,
  folder_output = "itsxpress_output",
  region = c("ITS1", "ITS2", "ALL"),
  taxa = "Fungi",
  paired_end = TRUE,
  single_file_output = FALSE,
  cluster_id = 1,
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
)
```

## Arguments

- path_to_fastq:

  (required) Path to the directory containing FASTQ files. See
  [`list_fastq_files()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/list_fastq_files.md)
  for help.

- folder_output:

  (default: `"itsxpress_output"`) Path to the output directory where
  trimmed FASTQ files and per-sample logs are written.

- region:

  (default: `"ITS1"`) The ITS region to keep: `"ITS1"`, `"ITS2"` or
  `"ALL"` (full ITS1 + 5.8S + ITS2).

- taxa:

  (default: `"Fungi"`) Taxonomic group sequenced (`--taxa`). One of
  `"Fungi"`, `"All"`, `"Alveolata"`, `"Bryophyta"`, `"Bacillariophyta"`,
  `"Amoebozoa"`, `"Euglenozoa"`, `"Chlorophyta"`, `"Rhodophyta"`,
  `"Phaeophyceae"`, `"Marchantiophyta"`, `"Metazoa"`, `"Oomycota"`,
  `"Haptophyceae"`, `"Raphidophyceae"`, `"Rhizaria"`, `"Synurophyceae"`,
  `"Tracheophyta"`, `"Eustigmatophyceae"`, or `"Parabasalia"`.

- paired_end:

  (logical, default `TRUE`) Whether the input FASTQ files are
  paired-end.

- single_file_output:

  (logical, default `FALSE`) Only used when `paired_end = TRUE`. If
  `FALSE` (default, required for downstream DADA2 paired-end denoising),
  `itsxpress` outputs unmerged R1 and R2 files via `--outfile` and
  `--outfile2`. If `TRUE`, `itsxpress` outputs a single merged FASTQ
  file per pair via `--outfile`.

- cluster_id:

  (numeric, default `1.0`) Percent identity in `[0.99, 1.0]` for
  clustering reads with VSEARCH (`--cluster_id`). Set to `1.0` (default)
  for exact dereplication when calling ASVs.

- reversed_primers:

  (logical, default `FALSE`) If `TRUE`, pass `--reversed_primers` for
  primers in reverse orientation (Taylor et al. 2016).

- trim_ccs:

  (logical, default `FALSE`) If `TRUE`, pass `--trim-ccs` for PacBio CCS
  long reads (only valid when `paired_end = FALSE`).

- nproc:

  (integer, default `1`) Number of processor threads (`--threads`).

- pattern:

  a pattern to filter files (passed on to list.files function).

- pattern_R1:

  a pattern to filter R1 files (default "*R1*")

- pattern_R2:

  a pattern to filter R2 files (default "*R2*")

- nb_files:

  the number of fastq files to list (default FALSE)

- cmd_is_run:

  (logical, default `TRUE`) Whether to execute the `itsxpress` commands.
  If `FALSE`, the function only builds and returns a named list of shell
  commands.

- return_file_path:

  (logical, default `FALSE`) If `TRUE` (when `cmd_is_run = TRUE`),
  return a list of output FASTQ file paths (`fnfs` and optionally
  `fnrs`), ready to be chained in a `targets` pipeline.

- itsxpresspath:

  (default:
  [`find_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_itsxpress.md))
  Path to the `itsxpress` executable.

- args_before_itsxpress:

  (String, default `""`) A one-line bash command run before `itsxpress`,
  e.g. the conda activation string returned by
  [`install_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsxpress.md).

- extra_itsxpress_args:

  (String, default `""`) Additional command-line flags passed verbatim
  to `itsxpress`.

## Value

If `cmd_is_run = FALSE`, a named list of `itsxpress` shell commands. If
`return_file_path = TRUE`, a named list of output file paths (`fnfs`,
and `fnrs` when `paired_end = TRUE` and `single_file_output = FALSE`).
Otherwise `NULL` invisibly, after running `itsxpress`.

## References

Rivers AR, Weber KC, Gardner TG, Liu S, Armstrong SD (2018). ITSxpress:
Software to rapidly trim internally transcribed spacer sequences with
quality scores for marker gene analysis. *F1000Research* 7: 1418.
[doi:10.12688/f1000research.15704.1](https://doi.org/10.12688/f1000research.15704.1)

Einarsson SV, Rivers AR (2024). ITSxpress version 2: software to rapidly
trim internal transcribed spacer sequences with quality scores for
amplicon sequencing. *Microbiology Spectrum* 12(4): e00601-24.
[doi:10.1128/spectrum.00601-24](https://doi.org/10.1128/spectrum.00601-24)

## See also

[`install_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsxpress.md),
[`is_itsxpress_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsxpress_installed.md),
[`find_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_itsxpress.md),
[`itsx_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsx_pq.md),
[`fastp()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/fastp.md),
[`cutadapt_remove_primers()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/cutadapt_remove_primers.md)

## Author

Adrien Taudière

## Examples

``` r
# Build paired-end ITSxpress commands without running them
itsxpress_pq(
  path_to_fastq = system.file("extdata", package = "MiscMetabar"),
  folder_output = file.path(tempdir(), "itsxpress_out"),
  region = "ITS2",
  cmd_is_run = FALSE
)
#> $`/tmp/RtmpwrUtcL/temp_libpath13913f132fe6a6/MiscMetabar/extdata/ex_R1_001.fastq.gz`
#> [1] "'itsxpress' --fastq '/tmp/RtmpwrUtcL/temp_libpath13913f132fe6a6/MiscMetabar/extdata/ex_R1_001.fastq.gz' --outfile '/tmp/Rtmp7j9BNw/itsxpress_out/ex_R1_001.fastq.gz' --region ITS2 --taxa Fungi --cluster_id 1 --threads 1 --log '/tmp/Rtmp7j9BNw/itsxpress_out/ex_R1_001.itsxpress.log' --fastq2 '/tmp/RtmpwrUtcL/temp_libpath13913f132fe6a6/MiscMetabar/extdata/ex_R2_001.fastq.gz' --outfile2 '/tmp/Rtmp7j9BNw/itsxpress_out/ex_R2_001.fastq.gz'"
#> 

# Run ITSxpress (the string returned by install_itsxpress() can be passed
# as `args_before_itsxpress` when it lives in a conda environment)
if (MiscMetabar::is_itsxpress_installed()) {
  itsxpress_pq(
    path_to_fastq = system.file("extdata", package = "MiscMetabar"),
    folder_output = file.path(tempdir(), "itsxpress_out"),
    region = "ITS1",
    taxa = "Fungi",
    nproc = 2
  )
}
```
