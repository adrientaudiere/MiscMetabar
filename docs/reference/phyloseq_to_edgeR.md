# Convert phyloseq OTU count data into DGEList for edgeR package

Convert phyloseq OTU count data into DGEList for edgeR package

## Usage

``` r
phyloseq_to_edgeR(physeq, group, method = "RLE", ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- group:

  (required) A character vector or factor giving the experimental
  group/condition for each sample/library. Alternatively, you may
  provide the name of a sample variable. This name should be among the
  output of `sample_variables(physeq)`, in which case
  `get_variable(physeq, group)` would return either a character vector
  or factor. This is passed on to
  [`DGEList`](https://rdrr.io/pkg/edgeR/man/DGEList.html), and you may
  find further details or examples in its documentation.

- method:

  The label of the edgeR-implemented normalization to use. See
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)
  for supported options and details. The default option is `"RLE"`,
  which is a scaling factor method proposed by Anders and Huber (2010).
  At time of writing, the
  [edgeR](https://rdrr.io/pkg/edgeR/man/edgeR-package.html) package
  supported the following options to the `method` argument:

  `c("TMM", "RLE", "upperquartile", "none")`.

- ...:

  Additional arguments passed on to
  [`DGEList`](https://rdrr.io/pkg/edgeR/man/DGEList.html)

## Value

A DGEList object. See
[`edgeR::estimateTagwiseDisp()`](https://rdrr.io/pkg/edgeR/man/estimateTagwiseDisp.html)
for more details.
