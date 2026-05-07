# Compute tSNE position of samples from a phyloseq object

Compute tSNE position of samples from a phyloseq object

## Usage

``` r
tsne_pq(physeq, method = "bray", dims = 2, theta = 0, perplexity = 30, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- method:

  A method to calculate distance using
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  function

- dims:

  (Int) Output dimensionality (default: 2)

- theta:

  (Numeric) Speed/accuracy trade-off (increase for less accuracy), set
  to 0.0 for exact TSNE (default: 0.0 see details in the man page of
  [`Rtsne::Rtsne`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html)).

- perplexity:

  (Numeric) Perplexity parameter (should not be bigger than 3 \*
  perplexity \< nrow(X) - 1, see details in the man page of
  [`Rtsne::Rtsne`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html))

- ...:

  Additional arguments passed on to
  [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html)

## Value

A list of element including the matrix Y containing the new
representations for the objects. See ?Rtsne::Rtsne() for more
information

## Examples

``` r
if (requireNamespace("Rtsne")) {
  res_tsne <- tsne_pq(data_fungi)
}
#> Taxa are now in rows.
```
