# Subset taxa using a taxa control or distribution based method

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

There is 3 main methods : discard taxa (i) using a control taxa (e.g.
truffle root tips), (ii) using a mixture models to detect bimodality in
pseudo-abundance distribution or (iii) using a minimum difference
threshold pseudo-abundance. Each cutoff is defined at the sample level.

## Usage

``` r
subset_taxa_tax_control(
  physeq,
  taxa_distri,
  method = "mean",
  min_diff_for_cutoff = NULL
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxa_distri:

  (required) a vector of length equal to the number of samples with the
  number of sequences per sample for the taxa control

- method:

  (default: "mean") a method to calculate the cut-off value. There are 6
  available methods:

  1.  `cutoff_seq`: discard taxa with less than the number of sequence
      than taxa control,

  2.  `cutoff_mixt`: using mixture models,

  3.  `cutoff_diff`: using a minimum difference threshold (need the
      argument min_diff_for_cutoff)

  4.  `min`: the minimum of the three firsts methods

  5.  `max`: the maximum of the three firsts methods

  6.  `mean`: the mean of the three firsts methods

- min_diff_for_cutoff:

  (int) argument for method `cutoff_diff`. Required if method is
  `cutoff_diff`, `min`, `max` or `mean`

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object.

## Author

Adrien Taudière

## Examples

``` r
subset_taxa_tax_control(data_fungi,
  as.numeric(data_fungi@otu_table[, 300]),
  min_diff_for_cutoff = 2
)
#> number of iterations= 6 
#> number of iterations= 19 
#> number of iterations= 7 
#> number of iterations= 19 
#> number of iterations= 5 
#> number of iterations= 14 
#> Error in stats::uniroot(f = f, lower = 1, upper = 1000) : 
#>   f.upper = f(upper) is NA
#> Warning: NAs introduced by coercion
#> number of iterations= 8 
#> number of iterations= 9 
#> number of iterations= 9 
#> number of iterations= 11 
#> number of iterations= 10 
#> number of iterations= 11 
#> number of iterations= 14 
#> number of iterations= 22 
#> number of iterations= 26 
#> number of iterations= 40 
#> number of iterations= 9 
#> number of iterations= 9 
#> number of iterations= 5 
#> number of iterations= 41 
#> number of iterations= 14 
#> number of iterations= 9 
#> number of iterations= 6 
#> number of iterations= 14 
#> number of iterations= 20 
#> number of iterations= 14 
#> number of iterations= 8 
#> number of iterations= 12 
#> number of iterations= 10 
#> number of iterations= 6 
#> number of iterations= 9 
#> number of iterations= 4 
#> number of iterations= 10 
#> number of iterations= 4 
#> number of iterations= 7 
#> number of iterations= 21 
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 9 
#> number of iterations= 8 
#> number of iterations= 19 
#> number of iterations= 10 
#> number of iterations= 8 
#> number of iterations= 21 
#> number of iterations= 6 
#> number of iterations= 15 
#> number of iterations= 17 
#> number of iterations= 3 
#> number of iterations= 26 
#> number of iterations= 11 
#> number of iterations= 11 
#> number of iterations= 15 
#> number of iterations= 5 
#> number of iterations= 23 
#> number of iterations= 9 
#> number of iterations= 23 
#> number of iterations= 12 
#> number of iterations= 7 
#> number of iterations= 11 
#> number of iterations= 25 
#> number of iterations= 24 
#> number of iterations= 12 
#> number of iterations= 4 
#> number of iterations= 15 
#> number of iterations= 10 
#> number of iterations= 6 
#> number of iterations= 7 
#> number of iterations= 12 
#> number of iterations= 10 
#> Error in stats::uniroot(f = f, lower = 1, upper = 1000) : 
#>   f.upper = f(upper) is NA
#> Warning: NAs introduced by coercion
#> number of iterations= 6 
#> number of iterations= 7 
#> number of iterations= 16 
#> number of iterations= 13 
#> number of iterations= 15 
#> number of iterations= 12 
#> number of iterations= 5 
#> number of iterations= 15 
#> number of iterations= 14 
#> number of iterations= 5 
#> number of iterations= 5 
#> Error in stats::uniroot(f = f, lower = 1, upper = 1000) : 
#>   f.upper = f(upper) is NA
#> Warning: NAs introduced by coercion
#> number of iterations= 5 
#> number of iterations= 10 
#> number of iterations= 17 
#> number of iterations= 3 
#> number of iterations= 14 
#> number of iterations= 14 
#> number of iterations= 9 
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 7 
#> number of iterations= 4 
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 6 
#> Warning: no non-missing arguments to min; returning Inf
#> number of iterations= 24 
#> number of iterations= 19 
#> number of iterations= 19 
#> number of iterations= 19 
#> number of iterations= 12 
#> number of iterations= 24 
#> number of iterations= 5 
#> number of iterations= 11 
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 6 
#> number of iterations= 22 
#> number of iterations= 11 
#> number of iterations= 8 
#> number of iterations= 13 
#> number of iterations= 42 
#> number of iterations= 4 
#> number of iterations= 7 
#> Error in stats::uniroot(f = f, lower = 1, upper = 1000) : 
#>   f.upper = f(upper) is NA
#> Warning: NAs introduced by coercion
#> number of iterations= 16 
#> number of iterations= 5 
#> number of iterations= 6 
#> number of iterations= 4 
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 4 
#> number of iterations= 14 
#> number of iterations= 17 
#> number of iterations= 22 
#> number of iterations= 6 
#> number of iterations= 6 
#> number of iterations= 3 
#> number of iterations= 28 
#> number of iterations= 10 
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to min; returning Inf
#> number of iterations= 3 
#> Warning: no non-missing arguments to min; returning Inf
#> number of iterations= 7 
#> number of iterations= 6 
#> number of iterations= 3 
#> number of iterations= 9 
#> number of iterations= 4 
#> number of iterations= 8 
#> number of iterations= 8 
#> number of iterations= 8 
#> Warning: no non-missing arguments to min; returning Inf
#> number of iterations= 8 
#> number of iterations= 24 
#> number of iterations= 37 
#> number of iterations= 13 
#> number of iterations= 3 
#> number of iterations= 34 
#> number of iterations= 3 
#> number of iterations= 5 
#> number of iterations= 6 
#> number of iterations= 6 
#> number of iterations= 12 
#> number of iterations= 5 
#> number of iterations= 7 
#> number of iterations= 9 
#> number of iterations= 12 
#> number of iterations= 14 
#> number of iterations= 10 
#> The filtering processes discard 61 taxa and 34991 sequences. Note that for  61 samples, all taxa were discarded. Please run clean_pq() to remove empty samples.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1359 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1359 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1359 reference sequences ]
```
