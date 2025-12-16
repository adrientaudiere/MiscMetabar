# Verify the validity of a phyloseq object

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Mostly for internal use in MiscMetabar functions.

## Usage

``` r
verify_pq(
  physeq,
  verbose = FALSE,
  min_nb_seq_sample = 500,
  min_nb_seq_taxa = 1
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- verbose:

  (logical, default FALSE) If TRUE, prompt some warnings.

- min_nb_seq_sample:

  (numeric) Only used if verbose = TRUE. Minimum number of sequences per
  samples to not show warning.

- min_nb_seq_taxa:

  (numeric) Only used if verbose = TRUE. Minimum number of sequences per
  taxa to not show warning.

## Value

Nothing if the phyloseq object is valid. An error in the other case.
Warnings if verbose = TRUE

## Author

Adrien Taudière
