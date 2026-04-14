# Variance Stabilizing Transformation of a phyloseq object (DESeq2)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Wrapper around
[`DESeq2::varianceStabilizingTransformation()`](https://rdrr.io/pkg/DESeq2/man/varianceStabilizingTransformation.html)
(Love, Huber & Anders 2014,
[doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
). Counts are incremented by 1 to handle zeros before VST is applied.

## Usage

``` r
vst_pq(physeq, blind = TRUE, fitType = "parametric")
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- blind:

  (logical, default `TRUE`) passed to DESeq2.

- fitType:

  (character, default `"parametric"`) passed to DESeq2.

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with a VST transformed `otu_table`.

## See also

[`DESeq2::varianceStabilizingTransformation()`](https://rdrr.io/pkg/DESeq2/man/varianceStabilizingTransformation.html)

## Author

Adrien Taudière

## Examples

``` r
data_f_vst <- vst_pq(data_fungi_mini)
#> converting counts to integer mode
#> -- note: fitType='parametric', but the dispersion trend was not well captured by the
#>    function: y = a/x + b, and a local regression fit was automatically substituted.
#>    specify fitType='local' or 'mean' to avoid this message next time.
```
