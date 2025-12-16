# Make upset plot for phyloseq object.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Alternative to venn plot.

## Usage

``` r
upset_pq(
  physeq,
  fact,
  taxa_fill = NULL,
  min_nb_seq = 0,
  na_remove = TRUE,
  numeric_fonction = sum,
  rarefy_after_merging = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor to cluster samples by modalities. Need
  to be in `physeq@sam_data`.

- taxa_fill:

  (default NULL) fill the ASV upset using a column in `tax_table` slot.

- min_nb_seq:

  minimum number of sequences by OTUs by samples to take into count this
  OTUs in this sample. For example, if min_nb_seq=2,each value of 2 or
  less in the OTU table will not count in the venn diagram

- na_remove:

  : if TRUE (the default), NA values in fact are removed if FALSE, NA
  values are set to "NA"

- numeric_fonction:

  (default : sum) the function for numeric vector useful only for
  complex plot (see examples)

- rarefy_after_merging:

  Rarefy each sample after merging by the modalities of `fact` parameter

- ...:

  Additional arguments passed on to the
  [`ComplexUpset::upset()`](https://krassowski.github.io/complex-upset/reference/upset.html)

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)2 plot

## See also

[`ggvenn_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggvenn_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("ComplexUpset") && packageVersion("ggplot2") < "4.0.0") {
  upset_pq(data_fungi_mini,
    fact = "Height", width_ratio = 0.2,
    taxa_fill = "Class"
  )
}
#> Loading required namespace: ComplexUpset
# \donttest{
if (requireNamespace("ComplexUpset") && packageVersion("ggplot2") < "4.0.0") {
  upset_pq(data_fungi_mini, fact = "Height", min_nb_seq = 1000)
  upset_pq(data_fungi_mini, fact = "Height", na_remove = FALSE)

  upset_pq(data_fungi_mini, fact = "Time", width_ratio = 0.2, rarefy_after_merging = TRUE)

  upset_pq(
    data_fungi_mini,
    fact = "Time",
    width_ratio = 0.2,
    annotations = list(
      "Sequences per ASV \n (log10)" = (
        ggplot(mapping = aes(y = log10(Abundance)))
        +
          geom_jitter(aes(
            color =
              Abundance
          ), na.rm = TRUE)
          +
          geom_violin(alpha = 0.5, na.rm = TRUE) +
          theme(legend.key.size = unit(0.2, "cm")) +
          theme(axis.text = element_text(size = 12))
      ),
      "ASV per phylum" = (
        ggplot(mapping = aes(fill = Phylum))
        +
          geom_bar() +
          ylab("ASV per phylum") +
          theme(legend.key.size = unit(0.2, "cm")) +
          theme(axis.text = element_text(size = 12))
      )
    )
  )

  upset_pq(
    data_fungi_mini,
    fact = "Time",
    width_ratio = 0.2,
    numeric_fonction = mean,
    annotations = list(
      "Sequences per ASV \n (log10)" = (
        ggplot(mapping = aes(y = log10(Abundance)))
        +
          geom_jitter(aes(
            color =
              Abundance
          ), na.rm = TRUE)
          +
          geom_violin(alpha = 0.5, na.rm = TRUE) +
          theme(legend.key.size = unit(0.2, "cm")) +
          theme(axis.text = element_text(size = 12))
      ),
      "ASV per phylum" = (
        ggplot(mapping = aes(fill = Phylum))
        +
          geom_bar() +
          ylab("ASV per phylum") +
          theme(legend.key.size = unit(0.2, "cm")) +
          theme(axis.text = element_text(size = 12))
      )
    )
  )

  upset_pq(
    subset_taxa(data_fungi_mini, Phylum == "Basidiomycota"),
    fact = "Time",
    width_ratio = 0.2,
    base_annotations = list(),
    annotations = list(
      "Sequences per ASV \n (log10)" = (
        ggplot(mapping = aes(y = log10(Abundance)))
        +
          geom_jitter(aes(
            color =
              Abundance
          ), na.rm = TRUE)
          +
          geom_violin(alpha = 0.5, na.rm = TRUE) +
          theme(legend.key.size = unit(0.2, "cm")) +
          theme(axis.text = element_text(size = 12))
      ),
      "ASV per phylum" = (
        ggplot(mapping = aes(fill = Class))
        +
          geom_bar() +
          ylab("ASV per Class") +
          theme(legend.key.size = unit(0.2, "cm")) +
          theme(axis.text = element_text(size = 12))
      )
    )
  )

  data_fungi2 <- data_fungi_mini
  data_fungi2@sam_data[["Time_0"]] <- data_fungi2@sam_data$Time == 0
  data_fungi2@sam_data[["Height__Time_0"]] <-
    paste0(data_fungi2@sam_data[["Height"]], "__", data_fungi2@sam_data[["Time_0"]])
  data_fungi2@sam_data[["Height__Time_0"]][grepl("NA", data_fungi2@sam_data[["Height__Time_0"]])] <-
    NA
  upset_pq(data_fungi2, fact = "Height__Time_0", width_ratio = 0.2, min_size = 2)
}
# }
```
