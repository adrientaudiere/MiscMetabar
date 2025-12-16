# Venn diagram of [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html) object

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Graphical representation of distribution of taxa across combined
modality of a factor.

## Usage

``` r
venn_pq(physeq, fact, min_nb_seq = 0, print_values = TRUE)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor to cluster samples by modalities. Need
  to be in `physeq@sam_data`.

- min_nb_seq:

  (default: 0) minimum number of sequences by OTUs by samples to take
  into count this OTUs in this sample. For example, if min_nb_seq=2,each
  value of 2 or less in the OTU table will be change into 0 for the
  analysis

- print_values:

  (logical) Print (or not) the table of number of OTUs for each
  combination. If print_values is TRUE the object is not a ggplot
  object. Please use print_values = FALSE if you want to add ggplot
  function (cf example).

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)2 plot
representing Venn diagram of modalities of the argument `factor`

## See also

[`venneuler`](https://rdrr.io/pkg/venneuler/man/venneuler.html)

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("venneuler")) {
  data("enterotype")
  venn_pq(enterotype, fact = "SeqTech")
}
#> Loading required namespace: venneuler

# \donttest{
if (requireNamespace("venneuler")) {
  venn_pq(enterotype, fact = "ClinicalStatus")
  venn_pq(enterotype, fact = "Nationality", print_values = FALSE)
  venn_pq(enterotype, fact = "ClinicalStatus", print_values = FALSE) +
    scale_fill_hue()
  venn_pq(enterotype, fact = "ClinicalStatus", print_values = FALSE) +
    scale_fill_hue()
}


# }
```
