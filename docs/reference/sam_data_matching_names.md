# Match sample names from sam_data and fastq files

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for targets bioinformatic pipeline.

## Usage

``` r
sam_data_matching_names(
  path_sam_data,
  sample_col_name,
  path_raw_seq,
  pattern_remove_sam_data = NULL,
  pattern_remove_fastq_files = NULL,
  verbose = TRUE,
  remove_undocumented_fastq_files = FALSE,
  prefix = NULL,
  ...
)
```

## Arguments

- path_sam_data:

  (Required) Path to sample data file.

- sample_col_name:

  (Required) The name of the column defining sample names in the sample
  data file.

- path_raw_seq:

  (Required) Path to the folder containing fastq files

- pattern_remove_sam_data:

  If not null, describe the pattern that will be deleted from sam_data
  samples names.

- pattern_remove_fastq_files:

  If not null, describe the pattern that will be deleted from fastq
  files names.

- verbose:

  (logical, default TRUE) If TRUE, print some additional messages.

- remove_undocumented_fastq_files:

  (logical, default FALSE) If set to TRUE fastq files not present in
  sam_data are removed from your folder. Keep a copy of those files
  somewhere before.

- prefix:

  Add a prefix to new samples names (ex. prefix = "samp")

- ...:

  Other parameters passed on to
  [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html)
  function.

## Value

A list of two objects :

- \$sam_names_matching is a tibble of corresponding samples names

- \$sam_data is a sample data files including only matching sample names

## Author

Adrien Taudière
