# Fastp: adapter detection, quality control and preprocessing of FASTQ files

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function is a wrapper for the
[fastp](https://github.com/OpenGene/fastp) software to perform quality
control and preprocessing of FASTQ files (paired-end or single-end).

To install fastp on Linux:

    wget http://opengene.org/fastp/fastp
    chmod a+x ./fastp
    sudo mv fastp /usr/local/bin/

For other OS or specific installation instructions, please refer to the
official fastp GitHub repository: <https://github.com/OpenGene/fastp>.

## Usage

``` r
fastp(
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
)
```

## Arguments

- path_to_fastq:

  (required) Path to the directory containing FASTQ files. See
  [`list_fastq_files()`](https://adrientaudiere.github.io/MiscMetabar/reference/list_fastq_files.md)
  for help.

- folder_output:

  Name of the output folder (default: "fastp_output").

- nproc:

  Number of threads to use (default: 1).

- pattern:

  a pattern to filter files (passed on to list.files function).

- pattern_R1:

  a pattern to filter R1 files (default "*R1*")

- pattern_R2:

  a pattern to filter R2 files (default "*R2*")

- nb_files:

  the number of fastq files to list (default FALSE)

- paired_end:

  Logical, whether files are paired-end (default: TRUE).

- cmd_is_run:

  Logical, whether to run the command (default: TRUE). If FALSE, the
  only effect of the function is to return a list of commands to run
  manually in a terminal.

- return_file_path:

  Logical, whether to return output file paths (default: FALSE). Useful
  in `targets` workflows.

- qualified_quality_phred:

  Minimum base quality (default: 20).

- average_qual:

  Minimum average read quality (default: 20).

- length_required:

  Minimum read length after trimming (default: 50).

- detect_adapter:

  Logical, enable automatic adapter detection for paired-end reads
  (default: TRUE).

- trim_poly_g:

  Logical, trim poly-G tails (default: TRUE).

- correction:

  Logical, correct mismatches in overlapped paired-end reads (default:
  TRUE).

- cut_front:

  Logical, trim low quality bases from 5' end (default: TRUE).

- cut_tail:

  Logical, trim low quality bases from 3' end (default: TRUE).

- cut_window_size:

  Window size for sliding quality filter (default: 4).

- cut_mean_quality:

  Minimum mean quality for trimming window (default: 20).

- extra_fastp_args:

  Additional fastp arguments as a character string (default: "").

## Value

If `cmd_is_run = FALSE`, a named list of fastp commands. If
`return_file_path = TRUE`, a list of output file paths. Otherwise `NULL`
invisibly, after running fastp.

## Details

This function is a wrapper for the fastp software, please cite fastp
([doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
) if you use this function.

## See also

[`is_fastp_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_fastp_installed.md),
[`list_fastq_files()`](https://adrientaudiere.github.io/MiscMetabar/reference/list_fastq_files.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# Paired-end example
fastp(
  path_to_fastq = system.file("extdata", package = "MiscMetabar"),
  folder_output = file.path(tempdir(), "qc_data"),
  nproc = 8,
  paired_end = TRUE
)

# Single-end example with extra options
fastp(
  path_to_fastq = "raw_data/",
  folder_output = file.path(tempdir(), "qc_data"),
  nproc = 8,
  paired_end = FALSE,
  extra_fastp_args = "--merge --n_base_limit 5 --dedup"
)

# Only build the commands, do not run them
fastp(
  path_to_fastq = system.file("extdata", package = "MiscMetabar"),
  cmd_is_run = FALSE
)
} # }
```
