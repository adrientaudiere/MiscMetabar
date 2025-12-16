# Add information to sample_data slot of a phyloseq-class object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Warning: The value nb_seq and nb_otu may be outdated if you transform
your phyloseq object, e.g. using the
[`subset_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/subset_taxa_pq.md)
function

## Usage

``` r
add_info_to_sam_data(
  physeq,
  df_info = NULL,
  add_nb_seq = TRUE,
  add_nb_otu = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- df_info:

  : A dataframe with rownames matching for sample names of the phyloseq
  object

- add_nb_seq:

  (Logical, default TRUE) Does we add a column nb_seq collecting the
  number of sequences per sample?

- add_nb_otu:

  (Logical, default TRUE) Does we add a column nb_otu collecting the
  number of OTUs per sample?

## Value

A phyloseq object with an updated sam_data slot

## Author

Adrien Taudière

## Examples

``` r
data_fungi <- add_info_to_sam_data(data_fungi)
boxplot(data_fungi@sam_data$nb_otu ~ data_fungi@sam_data$Time)


new_df <- data.frame(
  variable_1 = runif(n = nsamples(data_fungi), min = 1, max = 20),
  variable_2 = runif(n = nsamples(data_fungi), min = 1, max = 2)
)
rownames(new_df) <- sample_names(data_fungi)
data_fungi <- add_info_to_sam_data(data_fungi, new_df)
plot(data_fungi@sam_data$nb_otu ~ data_fungi@sam_data$variable_1)
```
