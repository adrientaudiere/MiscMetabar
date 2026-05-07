# Distribution of sequences across a factor for one taxon

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Focus on one taxon and one factor.

## Usage

``` r
distri_1_taxa(physeq, fact, taxa_name, digits = 2)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor in `physeq@sam_data` used to plot
  different lines

- taxa_name:

  (required) the name of the taxa

- digits:

  (default = 2) integer indicating the number of decimal places to be
  used (see [`?round`](https://rdrr.io/r/base/Round.html) for more
  information)

## Value

a dataframe with levels as rows and information as column :

- the number of sequences of the taxa (nb_seq)

- the number of samples of the taxa (nb_samp)

- the mean (mean_nb_seq) and standard deviation (sd_nb_seq) of the
  *nb_seq*

- the mean (mean_nb_seq_when_present) *nb_seq* excluding samples with
  zero

- the total number of samples (nb_total_samp)

- the proportion of samples with the taxa

## Author

Adrien Taudière

## Examples

``` r
distri_1_taxa(data_fungi, "Height", "ASV2")
#> Taxa are now in rows.
#>        nb_seq nb_samp mean_nb_seq sd_nb_seq mean_nb_seq_when_present
#> High    26876      28      655.51   2623.30                   959.86
#> Low     60925      29     1353.89   5075.12                  2100.86
#> Middle   4444      33       98.76    453.06                   134.67
#>        nb_total_samp prop_samp
#> High              41      0.68
#> Low               45      0.64
#> Middle            45      0.73
distri_1_taxa(data_fungi, "Time", "ASV81", digits = 1)
#> Taxa are now in rows.
#>    nb_seq nb_samp mean_nb_seq sd_nb_seq mean_nb_seq_when_present nb_total_samp
#> 0    7387       5       110.3     896.8                   1477.4            67
#> 5      24       5         0.6       1.9                      4.8            38
#> 10     10       3         0.5       1.4                      3.3            21
#> 15      0       0         0.0       0.0                      NaN            36
#>    prop_samp
#> 0        0.1
#> 5        0.1
#> 10       0.1
#> 15       0.0
```
