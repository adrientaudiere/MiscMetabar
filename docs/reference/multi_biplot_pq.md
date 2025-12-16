# Visualization of a collection of couples of samples for comparison

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This allow to plot all the possible
[`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md)
combination using one factor.

## Usage

``` r
multi_biplot_pq(physeq, split_by = NULL, pairs = NULL, na_remove = TRUE, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- split_by:

  (required if pairs is NULL) the name of the factor to make all
  combination of couples of values

- pairs:

  (required if split_by is NULL) the name of the factor in
  physeq@sam_data\` slot to make plot by pairs of samples. Each level
  must be present only two times. Note that if you set pairs, you also
  must set fact arguments to passed on to
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md).

- na_remove:

  (logical, default TRUE) if TRUE remove all the samples with NA in the
  `split_by` variable of the `physeq@sam_data` slot

- ...:

  Other parameters passed on to
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md)

## Value

a list of ggplot object

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
data_fungi_abun <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 10000)
#> Cleaning suppress 0 taxa (  ) and 15 sample(s) ( BE9-006-B_S27_MERGED.fastq.gz / C21-NV1-M_S64_MERGED.fastq.gz / DJ2-008-B_S87_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / DY5-004-M_S98_MERGED.fastq.gz / E9-009-B_S100_MERGED.fastq.gz / E9-009-H_S101_MERGED.fastq.gz / N22-001-B_S129_MERGED.fastq.gz / O20-X-B_S139_MERGED.fastq.gz / O21-007-M_S144_MERGED.fastq.gz / R28-008-H_S159_MERGED.fastq.gz / R28-008-M_S160_MERGED.fastq.gz / W26-001-M_S167_MERGED.fastq.gz / Y29-007-H_S182_MERGED.fastq.gz / Y29-007-M_S183_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1385
#> Number of kept ASV 35
#> Number of kept samples 170
p <- multi_biplot_pq(data_fungi_abun, "Height")
#> 54 were discarded due to NA in variables present in formula.
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
#> The two modalities differ greatly (more than x2) in their number of sequences (236336 vs 97106). You may be interested by the parameter rarefy_after_merging
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
#> Cleaning suppress 4 taxa and 0 samples.
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
lapply(p, print)



#> $`High-Low`

#> 
#> $`High-Middle`

#> 
#> $`Low-Middle`

#> 
# }
```
