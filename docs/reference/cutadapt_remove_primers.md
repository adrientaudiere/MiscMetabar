# Remove primers using [cutadapt](https://github.com/marcelm/cutadapt/)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

You need to install [Cutadapt](https://cutadapt.readthedocs.io/). See
also https://github.com/VascoElbrecht/JAMP/blob/master/JAMP/R/Cutadapt.R
for another call to cutadapt from R

## Usage

``` r
cutadapt_remove_primers(
  path_to_fastq,
  primer_fw = NULL,
  primer_rev = NULL,
  folder_output = "wo_primers",
  nproc = 1,
  pattern = "fastq.gz",
  pattern_R1 = "_R1",
  pattern_R2 = "_R2",
  nb_files = Inf,
  cmd_is_run = TRUE,
  return_file_path = FALSE,
  args_before_cutadapt =
    "source ~/miniconda3/etc/profile.d/conda.sh && conda activate cutadaptenv && "
)
```

## Arguments

- path_to_fastq:

  (Required) A path to a folder with fastq files. See
  [`list_fastq_files()`](https://adrientaudiere.github.io/MiscMetabar/reference/list_fastq_files.md)
  for help.

- primer_fw:

  (Required, String) The forward primer DNA sequence.

- primer_rev:

  (String) The reverse primer DNA sequence.

- folder_output:

  The path to a folder for output files

- nproc:

  (default 1) Set to number of cpus/processors to use for the clustering

- pattern:

  a pattern to filter files (passed on to list.files function).

- pattern_R1:

  a pattern to filter R1 files (default "*R1*")

- pattern_R2:

  a pattern to filter R2 files (default "*R2*")

- nb_files:

  the number of fastq files to list (default FALSE)

- cmd_is_run:

  (logical, default TRUE) Do the cutadapt command is run. If set to
  FALSE, the only effect of the function is to return a list of command
  to manually run in a terminal.

- return_file_path:

  (logical, default FALSE) If true, the function return the path of the
  output folder (param `folder_output`). Useful in targets workflow

- args_before_cutadapt:

  (String) A one line bash command to run before to run cutadapt. For
  examples, "source ~/miniconda3/etc/profile.d/conda.sh && conda
  activate cutadaptenv &&" allow to bypass the conda init which asks to
  restart the shell

## Value

a list of command or if `return_file_path` is TRUE, the path to the
output folder

## Details

This function is mainly a wrapper of the work of others. Please cite
cutadapt
([doi:10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200) ).

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
cutadapt_remove_primers(system.file("extdata", package = "MiscMetabar"),
  "TTC",
  "GAA",
  folder_output = tempdir()
)

cutadapt_remove_primers(
  system.file("extdata",
    package = "dada2"
  ),
  pattern_R1 = "F.fastq.gz",
  pattern_R2 = "R.fastq.gz",
  primer_fw = "TTC",
  primer_rev = "GAA",
  folder_output = tempdir()
)

cutadapt_remove_primers(
  system.file("extdata",
    package = "dada2"
  ),
  pattern_R1 = "F.fastq.gz",
  primer_fw = "TTC",
  folder_output = tempdir(),
  cmd_is_run = FALSE
)

unlink(tempdir(), recursive = TRUE)
} # }
```
