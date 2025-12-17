# Run Aldex on a phyloseq object

Run Aldex on a phyloseq object

## Usage

``` r
aldex_pq(physeq, bifactor = NULL, modalities = NULL, gamma = 0.5, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- bifactor:

  (required) The name of a column present in the `@sam_data` slot of the
  physeq object. Must be a character vector or a factor.

- modalities:

  (default NULL) A vector of modalities to keep in the analysis. If
  NULL, all modalities present in bifactor are kept. Note that only two
  modalities are allowed. @param gamma (default 0.5) The value of the
  Dirichlet Monte-Carlo sampling parameter. @param ... Additional
  arguments passed on to
  [`ALDEx2::aldex()`](https://rdrr.io/pkg/ALDEx2/man/aldex.html)

## Value

The result of
[`ALDEx2::aldex()`](https://rdrr.io/pkg/ALDEx2/man/aldex.html)

## Details

It is a wrapper of the
[`ALDEx2::aldex()`](https://rdrr.io/pkg/ALDEx2/man/aldex.html) function
with default gamma=0.5.
[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Author

Adrien Taudière

## Examples

``` r
res_aldex <- aldex_pq(data_fungi_mini,
  bifactor = "Height",
  modalities = c("Low", "High")
)
#> aldex.clr: generating Monte-Carlo instances and clr values
#> conditions vector supplied
#> operating in serial mode
#> aldex.scaleSim: adjusting samples to reflect scale uncertainty.
#> aldex.ttest: doing t-test
#> aldex.effect: calculating effect sizes
ALDEx2::aldex.plot(res_aldex, type = "volcano")
```
